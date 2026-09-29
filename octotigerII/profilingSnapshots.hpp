#pragma once

namespace octotigerII::profiling {
/// HPX startup callback, called on every locality before application setup.
/// Enabled by OCTOTIGERII_PROFILE_INTERVAL_SECONDS; inert without APEX.
void startSnapshots();
/// Console-locality call: stop and join snapshot callbacks on all localities.
void stopSnapshots();
}
