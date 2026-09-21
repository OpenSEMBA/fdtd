#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
#include <direct.h>
#include <windows.h>
#define mkdir_one(path) _mkdir(path)
#else
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>
#define mkdir_one(path) mkdir(path, 0777)
#endif

static int make_one_directory(const char *path) {
  if (mkdir_one(path) == 0 || errno == EEXIST) return 0;
  return errno;
}

int fdtd_create_directories(const char *path) {
  char *copy;
  size_t i;
  size_t length;
  int error;

  if (path == NULL || path[0] == '\0') return EINVAL;
  length = strlen(path);
  copy = malloc(length + 1);
  if (copy == NULL) return ENOMEM;
  memcpy(copy, path, length + 1);
  for (i = 1; i < length; ++i) {
    if (copy[i] == '/' || copy[i] == '\\') {
      char separator = copy[i];
      copy[i] = '\0';
      if (!(i == 2 && copy[1] == ':') && copy[0] != '\0') {
        error = make_one_directory(copy);
        if (error != 0) {
          free(copy);
          return error;
        }
      }
      copy[i] = separator;
    }
  }
  error = make_one_directory(copy);
  free(copy);
  return error;
}

int fdtd_delete_file(const char *path) {
  if (remove(path) == 0 || errno == ENOENT) return 0;
  return errno;
}

int fdtd_atomic_replace(const char *source, const char *target) {
  if (source == NULL || target == NULL || source[0] == '\0' || target[0] == '\0') {
    return EINVAL;
  }
#ifdef _WIN32
  if (MoveFileExA(source, target, MOVEFILE_REPLACE_EXISTING | MOVEFILE_WRITE_THROUGH)) {
    return 0;
  }
  return (int)GetLastError();
#else
  if (rename(source, target) == 0) return 0;
  return errno;
#endif
}

#ifdef _WIN32
int fdtd_remove_tree(const char *path) {
  char pattern[MAX_PATH];
  WIN32_FIND_DATAA entry;
  HANDLE handle;
  int error = 0;

  snprintf(pattern, sizeof(pattern), "%s\\*", path);
  handle = FindFirstFileA(pattern, &entry);
  if (handle == INVALID_HANDLE_VALUE) return GetLastError() == ERROR_FILE_NOT_FOUND ? 0 : (int)GetLastError();
  do {
    char child[MAX_PATH];
    if (strcmp(entry.cFileName, ".") == 0 || strcmp(entry.cFileName, "..") == 0) continue;
    snprintf(child, sizeof(child), "%s\\%s", path, entry.cFileName);
    if (entry.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) error = fdtd_remove_tree(child);
    else if (!DeleteFileA(child)) error = (int)GetLastError();
    if (error != 0) break;
  } while (FindNextFileA(handle, &entry));
  FindClose(handle);
  if (error == 0 && !RemoveDirectoryA(path)) error = (int)GetLastError();
  return error;
}
#else
int fdtd_remove_tree(const char *path) {
  DIR *directory;
  struct dirent *entry;
  int error = 0;

  directory = opendir(path);
  if (directory == NULL) return errno == ENOENT ? 0 : errno;
  while ((entry = readdir(directory)) != NULL) {
    char *child;
    size_t length;
    struct stat status;
    if (strcmp(entry->d_name, ".") == 0 || strcmp(entry->d_name, "..") == 0) continue;
    length = strlen(path) + strlen(entry->d_name) + 2;
    child = malloc(length);
    if (child == NULL) {
      error = ENOMEM;
      break;
    }
    snprintf(child, length, "%s/%s", path, entry->d_name);
    if (lstat(child, &status) != 0) error = errno;
    else if (S_ISDIR(status.st_mode)) error = fdtd_remove_tree(child);
    else if (unlink(child) != 0) error = errno;
    free(child);
    if (error != 0) break;
  }
  closedir(directory);
  if (error == 0 && rmdir(path) != 0) error = errno;
  return error;
}
#endif
