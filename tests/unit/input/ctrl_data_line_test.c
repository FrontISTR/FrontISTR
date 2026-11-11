/*****************************************************************************
 * Copyright (c) 2026 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#include "fstr_ctrl_util.h"

static int check_line(const char *filename, const char *text, int inp,
                      int expected_line) {
  FILE *fp;
  fstr_ctrl_data *ctrl;
  int ids[2] = {0, 0};
  double values[2] = {0.0, 0.0};
  char message[1024], expected[128];
  int failed = 0;

  fp = fopen(filename, "w");
  if (!fp) return 1;
  fputs(text, fp);
  if (fclose(fp)) return 1;
  ctrl = inp ? c_fstr_ctrl_open_inp(filename) : c_fstr_ctrl_open(filename);
  if (!ctrl) return 1;

  snprintf(expected, sizeof(expected), "line:%d, column:2", expected_line);
  if (c_fstr_ctrl_get_data_array_ex(ctrl, "IR", ids, values) == 0) failed = 1;
  c_fstr_ctrl_get_err_msg(message, sizeof(message));
  if (!strstr(message, expected)) {
    fprintf(stderr, "%s: expected %s, got %s", filename, expected, message);
    failed = 1;
  }
  if (c_fstr_ctrl_get_data_ex(ctrl, 2, "IR", &ids[0], &values[0]) == 0)
    failed = 1;
  c_fstr_ctrl_get_err_msg(message, sizeof(message));
  if (!strstr(message, expected)) {
    fprintf(stderr, "%s (single row): expected %s, got %s",
            filename, expected, message);
    failed = 1;
  }
  if (c_fstr_ctrl_get_data_ex(ctrl, 3, "IR", &ids[0], &values[0]) == 0)
    failed = 1;
  c_fstr_ctrl_get_err_msg(message, sizeof(message));
  if (!strstr(message, "line:2,")) {
    fprintf(stderr, "%s (missing row): expected header line 2, got %s",
            filename, message);
    failed = 1;
  }
  c_fstr_ctrl_close(ctrl);
  remove(filename);
  return failed;
}

int main(void) {
  int failed = 0;
  failed += check_line("ctrl_line_cnt.ctl",
      "!! control input\n!CLOAD\n!! comment\n# comment\n\n"
      "1,1.0\n!! comment\n2,bad\n", 0, 8);
  failed += check_line("ctrl_line_inp.ctl",
      "** control input\n*CLOAD,\n AMPLITUDE=RAMP\n** comment\n\n"
      "1,1.0\n** comment\n\n2,bad\n", 1, 9);
  if (failed) return EXIT_FAILURE;
  puts("cnt and inp errors use physical data lines; missing rows use the header");
  return EXIT_SUCCESS;
}
