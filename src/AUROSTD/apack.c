// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
// This is a minimal tool to create xz compressed tar archives to prepare data
// that fot embedding in src/AUROSTD/aurostd_data.cpp
// note: should only be utilized by the build process
// usage: apack ARCHIVE_NAME [FILES_TO_PACK ...]
// author: hagen.eckert@duke.edu

#include <archive.h>
#include <archive_entry.h>

#include <fcntl.h>
#include <stdio.h>
#include <unistd.h>

static char buff[16384];

int main(int argc, const char **argv) {
  struct archive *a;
  struct archive *disk;
  struct archive_entry *entry;
  int r;
  int fd;
  ssize_t len;

  // initialized the new archive
  a = archive_write_new();
  archive_write_set_format_pax_restricted(a);
  archive_write_add_filter_xz(a);
  archive_write_set_options(a, "compression-level=5,threads=1");

  // open the archive to write
  r = archive_write_open_filename(a, argv[1]);
  if (r != ARCHIVE_OK) {
    printf("%s\n", archive_error_string(a));
    return (1);
  }

  // iterate over files to pack
  for (int i = 2; i < argc; ++i) {
    disk = archive_read_disk_new();
    archive_read_disk_set_standard_lookup(disk);
    r = archive_read_disk_open(disk, argv[i]);
    if (r != ARCHIVE_OK) {
      printf("%s\n", archive_error_string(a));
      return (1);
    }
    for (;;) {
      entry = archive_entry_new();
      r = archive_read_next_header2(disk, entry);
      if (r == ARCHIVE_EOF) {
        break;
      }
      if (r != ARCHIVE_OK) {
        printf("%s\n", archive_error_string(disk));
        return (1);
      }
      archive_read_disk_descend(disk);
      r = archive_write_header(a, entry);
      if (r < ARCHIVE_OK) {
        printf("%s\n", archive_error_string(a));
      }
      if (r == ARCHIVE_FATAL) {
        printf("FATAL error - unable to write entry to archive.");
        return (1);
      }
      if (r > ARCHIVE_FAILED) {
        fd = open(archive_entry_sourcepath(entry), O_RDONLY);
        len = read(fd, buff, sizeof(buff));
        while (len > 0) {
          archive_write_data(a, buff, len);
          len = read(fd, buff, sizeof(buff));
        }
        close(fd);
      }
      archive_entry_free(entry);
    }
    archive_write_free(disk);
  }
  archive_write_close(a);
  archive_write_free(a);
  return (0);
}
