"""Run a test executable on independent HPX TCP localities on this host.

Usage: python3 tests/distributed.py build/tests/storageChecks [application args]
This checks distributed correctness, not cluster performance.
"""
import os
from pathlib import Path
import socket
import subprocess
import sys
import tempfile


def run(executable, arguments, localities=2):
    if localities < 1:
        raise ValueError("At least one HPX locality is required")
    # Keep port reservations until all numbers have been selected.
    sockets = [socket.socket() for _ in range(localities)]
    for sock in sockets:
        sock.bind(("127.0.0.1", 0))
    ports = [sock.getsockname()[1] for sock in sockets]
    for sock in sockets:
        sock.close()
    executable = str(Path(executable).resolve())
    common = [f"--hpx:localities={localities}", "--hpx:threads=" + os.environ.get("OCTOTIGERII_TEST_THREADS", "2"), "--hpx:bind=none",
              f"--hpx:agas=127.0.0.1:{ports[0]}",
              "--hpx:ini=hpx.parcel.message_handlers=" + os.environ.get("OCTOTIGERII_TEST_COALESCING", "1"),
              "--hpx:ini=hpx.parcel.zero_copy_serialization_threshold=" + os.environ.get("OCTOTIGERII_TEST_CHUNK_THRESHOLD", "1")]
    with tempfile.TemporaryDirectory(prefix="octotigerII-distributed-") as folder:
        logs = [open(Path(folder) / f"locality-{i}.log", "w+") for i in range(localities)]
        processes = []
        try:
            for i in range(localities):
                command = [executable, *arguments, *common,
                           f"--hpx:hpx=127.0.0.1:{ports[i]}"]
                if i:
                    command.append("--hpx:worker")
                processes.append(subprocess.Popen(command, stdout=logs[i], stderr=subprocess.STDOUT))
            statuses = [process.wait(timeout=300) for process in processes]
            if any(statuses):
                raise RuntimeError(f"Distributed execution failed: {statuses}")
        finally:
            for process in processes:
                if process.poll() is None:
                    process.kill()
                    process.wait()
            for i, log in enumerate(logs):
                log.seek(0)
                print(f"Locality {i}:\n{log.read()}", end="")
                log.close()


if __name__ == "__main__":
    count = int(os.environ.get("OCTOTIGERII_TEST_LOCALITIES", "2"))
    run(sys.argv[1], sys.argv[2:], count)
