/* Operating-system memory probe used by the startup memory check.
 *
 * All functions return a byte count, or -1 when the value cannot be
 * determined on the current platform.  Keeping this in C avoids spreading
 * POSIX/Win32 headers through the Fortran sources.
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
#include <windows.h>
#include <psapi.h>
#else
#include <unistd.h>
#include <sys/resource.h>
#endif

#define FDTD_UNKNOWN ((int64_t)-1)

#if defined(__linux__)

static int64_t read_meminfo_kb(const char *key)
{
  FILE *file;
  char line[256];
  char label[64];
  unsigned long long value;
  size_t key_length;

  file = fopen("/proc/meminfo", "r");
  if (file == NULL) return FDTD_UNKNOWN;

  key_length = strlen(key);
  while (fgets(line, sizeof(line), file) != NULL) {
    if (strncmp(line, key, key_length) != 0) continue;
    if (sscanf(line, "%63[^:]: %llu", label, &value) != 2) continue;
    fclose(file);
    return (int64_t)value * 1024;
  }

  fclose(file);
  return FDTD_UNKNOWN;
}

static int64_t read_first_line_int64(const char *path)
{
  FILE *file;
  char line[128];
  char *end;
  long long value;

  file = fopen(path, "r");
  if (file == NULL) return FDTD_UNKNOWN;
  if (fgets(line, sizeof(line), file) == NULL) {
    fclose(file);
    return FDTD_UNKNOWN;
  }
  fclose(file);

  if (strncmp(line, "max", 3) == 0) return FDTD_UNKNOWN;

  value = strtoll(line, &end, 10);
  if (end == line) return FDTD_UNKNOWN;
  return (int64_t)value;
}

static int is_unlimited_cgroup_value(int64_t value, int64_t total)
{
  if (value < 0) return 1;
  if (value > (INT64_C(1) << 60)) return 1;
  if (total > 0 && value > total * 4) return 1;
  return 0;
}

static int64_t cgroup_headroom_for_path(const char *base, const char *path)
{
  char limit_path[1024];
  char usage_path[1024];
  int64_t limit;
  int64_t usage;

  snprintf(limit_path, sizeof(limit_path), "%s%s/memory.max", base, path);
  snprintf(usage_path, sizeof(usage_path), "%s%s/memory.current", base, path);
  limit = read_first_line_int64(limit_path);
  usage = read_first_line_int64(usage_path);
  if (limit != FDTD_UNKNOWN && usage != FDTD_UNKNOWN && usage >= 0) {
    if (!is_unlimited_cgroup_value(limit, 0)) return limit - usage;
  }

  snprintf(limit_path, sizeof(limit_path), "%s/memory%s/memory.limit_in_bytes", base, path);
  snprintf(usage_path, sizeof(usage_path), "%s/memory%s/memory.usage_in_bytes", base, path);
  limit = read_first_line_int64(limit_path);
  usage = read_first_line_int64(usage_path);
  if (limit != FDTD_UNKNOWN && usage != FDTD_UNKNOWN && usage >= 0) {
    if (!is_unlimited_cgroup_value(limit, 0)) return limit - usage;
  }

  return FDTD_UNKNOWN;
}

static int64_t cgroup_available_bytes(void)
{
  FILE *file;
  char line[512];
  char path[768];
  int64_t headroom;
  int64_t root_headroom;

  /* cgroup v2: "0::/path". cgroup v1: "...:memory:/path". */
  path[0] = '\0';
  file = fopen("/proc/self/cgroup", "r");
  if (file != NULL) {
    while (fgets(line, sizeof(line), file) != NULL) {
      char *separator = strchr(line, ':');
      if (separator == NULL) continue;
      if (strncmp(line, "0::", 3) == 0) {
        snprintf(path, sizeof(path), "%s", line + 3);
        break;
      }
      separator = strchr(separator + 1, ':');
      if (separator == NULL) continue;
      if (strncmp(line, "1:memory:", 9) == 0) {
        snprintf(path, sizeof(path), "%s", separator + 1);
        break;
      }
    }
    fclose(file);
  }

  {
    char *newline = strchr(path, '\n');
    if (newline != NULL) *newline = '\0';
  }

  headroom = FDTD_UNKNOWN;
  if (path[0] != '\0') {
    headroom = cgroup_headroom_for_path("/sys/fs/cgroup", path);
  }
  if (headroom == FDTD_UNKNOWN) {
    headroom = cgroup_headroom_for_path("/sys/fs/cgroup", "");
  }

  root_headroom = cgroup_headroom_for_path("/sys/fs/cgroup", "");
  if (root_headroom != FDTD_UNKNOWN &&
      (headroom == FDTD_UNKNOWN || root_headroom < headroom)) {
    headroom = root_headroom;
  }

  if (headroom != FDTD_UNKNOWN && headroom < 0) headroom = 0;
  return headroom;
}

static int64_t read_status_kb(const char *key)
{
  FILE *file;
  char line[256];
  unsigned long long value;

  file = fopen("/proc/self/status", "r");
  if (file == NULL) return FDTD_UNKNOWN;
  while (fgets(line, sizeof(line), file) != NULL) {
    if (strncmp(line, key, strlen(key)) != 0) continue;
    if (sscanf(line, "%*[^:]: %llu", &value) != 1) continue;
    fclose(file);
    return (int64_t)value * 1024;
  }
  fclose(file);
  return FDTD_UNKNOWN;
}

int64_t fdtd_memory_total_bytes(void)
{
  return read_meminfo_kb("MemTotal");
}

int64_t fdtd_memory_available_bytes(void)
{
  int64_t available;
  int64_t cgroup;

  available = read_meminfo_kb("MemAvailable");
  if (available == FDTD_UNKNOWN) {
    int64_t free_bytes = read_meminfo_kb("MemFree");
    int64_t buffers = read_meminfo_kb("Buffers");
    int64_t cached = read_meminfo_kb("Cached");
    if (free_bytes != FDTD_UNKNOWN && buffers != FDTD_UNKNOWN &&
        cached != FDTD_UNKNOWN) {
      available = free_bytes + buffers + cached;
    }
  }

  cgroup = cgroup_available_bytes();
  if (available == FDTD_UNKNOWN) return cgroup;
  if (cgroup != FDTD_UNKNOWN && cgroup < available) return cgroup;
  return available;
}

int64_t fdtd_process_resident_bytes(void)
{
  return read_status_kb("VmRSS:");
}

int64_t fdtd_process_virtual_bytes(void)
{
  int64_t virtual_bytes = read_status_kb("VmSize:");
  int64_t data_bytes = read_status_kb("VmData:");
  if (virtual_bytes == FDTD_UNKNOWN) return data_bytes;
  if (data_bytes == FDTD_UNKNOWN) return virtual_bytes;
  return virtual_bytes > data_bytes ? virtual_bytes : data_bytes;
}

int64_t fdtd_process_memory_limit_bytes(void)
{
  struct rlimit address_limit;
  struct rlimit data_limit;
  int64_t limit = FDTD_UNKNOWN;

  if (getrlimit(RLIMIT_AS, &address_limit) == 0 &&
      address_limit.rlim_cur != RLIM_INFINITY) {
    limit = (int64_t)address_limit.rlim_cur;
  }
  if (getrlimit(RLIMIT_DATA, &data_limit) == 0 &&
      data_limit.rlim_cur != RLIM_INFINITY) {
    int64_t data_bytes = (int64_t)data_limit.rlim_cur;
    if (limit == FDTD_UNKNOWN || data_bytes < limit) limit = data_bytes;
  }
  if (limit == FDTD_UNKNOWN) return 0;
  return limit;
}

int fdtd_host_name(char *buffer, size_t size)
{
  if (buffer == NULL || size == 0) return -1;
  if (gethostname(buffer, size) != 0) {
    buffer[0] = '\0';
    return -1;
  }
  buffer[size - 1] = '\0';
  return 0;
}

#elif defined(_WIN32)

int64_t fdtd_memory_total_bytes(void)
{
  MEMORYSTATUSEX status;
  status.dwLength = sizeof(status);
  if (!GlobalMemoryStatusEx(&status)) return FDTD_UNKNOWN;
  return (int64_t)status.ullTotalPhys;
}

int64_t fdtd_memory_available_bytes(void)
{
  MEMORYSTATUSEX status;
  int64_t available;
  int64_t commit_available;

  status.dwLength = sizeof(status);
  if (!GlobalMemoryStatusEx(&status)) return FDTD_UNKNOWN;

  available = (int64_t)status.ullAvailPhys;
  commit_available = (int64_t)status.ullAvailPageFile;
  if (commit_available > 0 && commit_available < available) {
    available = commit_available;
  }
  return available;
}

int64_t fdtd_process_resident_bytes(void)
{
  PROCESS_MEMORY_COUNTERS counters;
  if (!GetProcessMemoryInfo(GetCurrentProcess(), &counters, sizeof(counters))) {
    return FDTD_UNKNOWN;
  }
  return (int64_t)counters.WorkingSetSize;
}

int64_t fdtd_process_virtual_bytes(void)
{
  PROCESS_MEMORY_COUNTERS counters;
  if (!GetProcessMemoryInfo(GetCurrentProcess(), &counters, sizeof(counters))) {
    return FDTD_UNKNOWN;
  }
  return (int64_t)counters.PagefileUsage;
}

int64_t fdtd_process_memory_limit_bytes(void)
{
  BOOL in_job = FALSE;
  JOBOBJECT_EXTENDED_LIMIT_INFORMATION information;

  if (!IsProcessInJob(GetCurrentProcess(), NULL, &in_job) || !in_job) {
    return 0;
  }
  if (!QueryInformationJobObject(NULL, JobObjectExtendedLimitInformation,
                                 &information, sizeof(information), NULL)) {
    return 0;
  }
  if ((information.BasicLimitInformation.LimitFlags &
       JOB_OBJECT_LIMIT_PROCESS_MEMORY) == 0) {
    return 0;
  }
  return (int64_t)information.ProcessMemoryLimit;
}

int fdtd_host_name(char *buffer, size_t size)
{
  DWORD length;
  if (buffer == NULL || size == 0) return -1;
  length = (DWORD)size;
  if (!GetComputerNameA(buffer, &length)) {
    buffer[0] = '\0';
    return -1;
  }
  buffer[size - 1] = '\0';
  return 0;
}

#else

int64_t fdtd_memory_total_bytes(void) { return FDTD_UNKNOWN; }
int64_t fdtd_memory_available_bytes(void) { return FDTD_UNKNOWN; }
int64_t fdtd_process_resident_bytes(void) { return FDTD_UNKNOWN; }
int64_t fdtd_process_virtual_bytes(void) { return FDTD_UNKNOWN; }
int64_t fdtd_process_memory_limit_bytes(void) { return 0; }

int fdtd_host_name(char *buffer, size_t size)
{
  if (buffer == NULL || size == 0) return -1;
  buffer[0] = '\0';
  return -1;
}

#endif
