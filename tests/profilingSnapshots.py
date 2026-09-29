"""Verify live, cumulative APEX snapshots, abrupt death, and normal shutdown.

Uses real OctoII SCF work, with a short wall interval only for this test.
"""
import os
from pathlib import Path
import re
import socket
import subprocess
import sys
import tempfile
import time


def launch(executable, folder, count, interval, cells):
    reservations = [socket.socket() for _ in range(count)]
    for sock in reservations:
        sock.bind(("127.0.0.1", 0))
    ports = [sock.getsockname()[1] for sock in reservations]
    for sock in reservations:
        sock.close()
    processes, logs = [], []
    for rank in range(count):
        root = folder / f"rank-{rank}"
        root.mkdir()
        env = dict(os.environ, APEX_OUTPUT_FILE_PATH=str(root),
                   APEX_PROFILE_OUTPUT="1", APEX_SCREEN_OUTPUT="1", APEX_CSV_OUTPUT="1",
                   APEX_FINAL_OUTPUT_ONLY="0", APEX_DISABLE="0",
                   OCTOTIGERII_PROFILE_INTERVAL_SECONDS=str(interval))
        log = open(folder / f"rank-{rank}.log", "w+")
        args = [executable, "--problem.name=binary-scf", f"--scf.cells={cells}",
                "--scf.referenceWidth=2", "--scf.commonPolytropicK=on", "--scf.massRatio=.7",
                "--scf.virialTolerance=1", "--mesh.cells=4", "--mesh.level=1",
                "--runtime.stopTime=0", "--output.enabled=off", "--verification.analytic=off",
                f"--output.directory={folder / 'simulation'}", "--hpx:threads=2", "--hpx:bind=none",
                f"--hpx:localities={count}", f"--hpx:agas=127.0.0.1:{ports[0]}",
                f"--hpx:hpx=127.0.0.1:{ports[rank]}"]
        if rank:
            args.append("--hpx:worker")
        processes.append(subprocess.Popen(args, env=env, stdout=log, stderr=subprocess.STDOUT))
        logs.append(log)
    return processes, logs


def snapshots(folder, rank):
    return sorted((folder / f"rank-{rank}" / "snapshots" / f"locality-{rank}").glob("*.profile"))


def calls(text):
    match = re.search(r'^"scf.gravity_fft"\s+(\d+)', text, re.M)
    return int(match[1]) if match else 0


def cleanup(processes, logs):
    for process in processes:
        if process.poll() is None:
            process.kill()
        process.wait(timeout=30)
    for log in logs:
        log.close()


def test(executable, count):
    with tempfile.TemporaryDirectory(prefix="octo-apex-snapshots-") as temp:
        folder = Path(temp)
        live = folder / "killed"
        live.mkdir()
        processes, logs = launch(executable, live, count, .1, 64)
        try:
            deadline = time.monotonic() + 90
            while True:
                assert all(p.poll() is None for p in processes), "exited before live snapshot checks"
                try:
                    saved = [snapshots(live, rank) for rank in range(count)]
                    contents = [[f.read_text() for f in files] for files in saved]
                    ready = all(len(rows) == 2 and calls(rows[0]) >= 2 and
                                calls(rows[1]) >= calls(rows[0]) for rows in contents)
                except FileNotFoundError:  # retention raced this observer; retry
                    ready = False
                if ready:
                    break
                assert time.monotonic() < deadline, "no cumulative snapshots on every locality"
                time.sleep(.05)
            # No finalize handlers can run after SIGKILL. Completed files must
            # remain usable, even if the kill interrupts another profile write.
            for process in processes:
                process.kill()
            for process in processes:
                process.wait(timeout=30)
            for rank in range(count):
                files = snapshots(live, rank)
                assert 2 <= len(files) <= 3  # a kill can precede retention cleanup
                assert all("templated_functions_MULTI_TIME" in f.read_text() for f in files)
                assert max(calls(f.read_text()) for f in files) >= 2
            print(f"Live snapshots and SIGKILL survival passed on {count} localities")
        except Exception:
            for log in logs:
                log.flush()
                print(Path(log.name).read_text())
            raise
        finally:
            cleanup(processes, logs)
        for interval in (.1, 0):
            normal = folder / f"normal-{interval}"
            normal.mkdir()
            processes, logs = launch(executable, normal, count, interval, 32)
            try:
                for process in processes:
                    assert process.wait(timeout=90) == 0, "normal shutdown failed"
                for rank in range(count):
                    root = normal / f"rank-{rank}"
                    final = (root / f"profile.{rank}.0.0").read_text()
                    assert calls(final) > 0, "final APEX output settings not restored"
                    files = snapshots(normal, rank)
                    if interval:
                        assert len(files) == 2
                        assert calls(final) >= calls(files[-1].read_text()), "dump reset cumulative counters"
                    else:
                        assert not files, "disabled snapshots produced files"
                    assert (root / "apex_profiles.csv").exists(), "final CSV output missing"
                print(f"Normal exit and interval={interval} passed on {count} localities")
            finally:
                cleanup(processes, logs)


if __name__ == "__main__":
    test(str(Path(sys.argv[1]).resolve()), int(sys.argv[2]))
