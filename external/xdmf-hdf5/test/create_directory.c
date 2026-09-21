#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#ifdef _WIN32
#include <direct.h>
#define mkdir(path, mode) _mkdir(path)
#endif

static int create_one_directory(const char *path) {
  struct stat info;

  if (mkdir(path, 0777) != 0 && errno != EEXIST) return errno;
  if (stat(path, &info) != 0) return errno;
  return S_ISDIR(info.st_mode) ? 0 : ENOTDIR;
}

int xdmf_create_directory(const char *directory) {
  char *path;
  size_t i, length;
  int error;

  if (directory == NULL || directory[0] == '\0') return EINVAL;
  length = strlen(directory);
  path = malloc(length + 1);
  if (path == NULL) return ENOMEM;
  memcpy(path, directory, length + 1);

  for (i = 1; i < length; ++i) {
    if (path[i] != '/' && path[i] != '\\') continue;
#ifdef _WIN32
    if (i == 2 && path[1] == ':') continue;
#endif
    path[i] = '\0';
    error = create_one_directory(path);
    path[i] = directory[i];
    if (error != 0) {
      free(path);
      return error;
    }
  }

  error = create_one_directory(path);
  free(path);
  return error;
}
