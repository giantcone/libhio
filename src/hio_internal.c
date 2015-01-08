/* -*- Mode: C; c-basic-offset:2 ; indent-tabs-mode:nil -*- */
/*
 * Copyright (c) 2014      Los Alamos National Security, LLC.  All rights
 *                         reserved. 
 * $COPYRIGHT$
 * 
 * Additional copyrights may follow
 * 
 * $HEADER$
 */

#include "hio_types.h"

#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

#include <errno.h>
#include <pthread.h>
#include <unistd.h>
#include <string.h>

typedef struct hio_error_stack_item_t {
  struct hio_error_stack_item_t *next;
  hio_object_t                   object;
  int                            hrc;
  char                          *error_string;
} hio_error_stack_item_t;

static hio_error_stack_item_t *hio_error_stack_head = NULL;
static pthread_mutex_t hio_error_stack_mutex = PTHREAD_MUTEX_INITIALIZER;

/**
 * @file Internal hio functions
 */

void hioi_log (hio_context_t context, int level, char *format, ...) {
  if (context->context_verbose >= level) {
    time_t current_time;
    char time_buf[26];
    va_list vargs;

    current_time = time (NULL);
    ctime_r (&current_time, time_buf);
    time_buf[strlen(time_buf) - 1] = '\0';

    va_start (vargs, format);
    fprintf (stderr, "[hio:%d] (context: %s) %s: ", level, context->context_object.identifier, time_buf);
    vfprintf (stderr, format, vargs);
    fprintf (stderr, "\n");
    va_end (vargs);
  }
}

int hioi_err_errno (int err) {
  switch (err) {
  case EPERM:
  case EACCES:
    return HIO_ERR_PERM;
  case ENOMEM:
    return HIO_ERR_OUT_OF_RESOURCE;
  case ENOENT:
    return HIO_ERR_NOT_FOUND;
  case EIO:
    return HIO_ERR_IO_PERMANENT;
  default:
    return HIO_ERROR;
  }
}

void hio_err_push (int hrc, hio_context_t context, hio_object_t object, char *format, ...) {
  hio_error_stack_item_t *new_item;
  va_list vargs;
  int rc;

  new_item = calloc (1, sizeof (hio_error_stack_item_t));
  if (NULL == new_item) {
    /* not much can be done here. we are just plain OOM. */
    return;
  }

  va_start (vargs, format);

  rc = vasprintf (&new_item->error_string, format, vargs);

  va_end (vargs);

  if (0 >= rc) {
    /* couldn't allocate error string */
    free (new_item);
    return;
  }

  if (context) {
    hioi_log (context, HIO_VERBOSE_ERROR, "%s", new_item->error_string);
  }

  new_item->hrc = hrc;

  /* push the error message onto the stack */
  if (NULL == context) {
    pthread_mutex_lock (&hio_error_stack_mutex);
    new_item->next = hio_error_stack_head;
    hio_error_stack_head = new_item;
    pthread_mutex_unlock (&hio_error_stack_mutex);
  } else {
    pthread_mutex_lock (&context->context_lock);
    new_item->next = (hio_error_stack_item_t *) context->context_error_stack;
    context->context_error_stack = (void *) new_item;
    pthread_mutex_unlock (&context->context_lock);
  }
}

#if HIO_USE_MPI
void hio_err_push_mpi (int mpirc, hio_context_t context, hio_object_t object, char *format, ...) {
  hio_error_stack_item_t *new_item;
  char mpi_error[MPI_MAX_ERROR_STRING] = "Unknown error";
  int resultlen = MPI_MAX_ERROR_STRING;
  va_list vargs;
  char *temp;
  int rc;

  va_start (vargs, format);

  rc = vasprintf (&temp, format, vargs);

  va_end (vargs);

  if (0 >= rc) {
    /* couldn't allocate error string */
    free (new_item);
    return;
  }

  /* ignore the error code for this */
  (void) MPI_Error_string (mpirc, mpi_error, &resultlen);

  new_item = calloc (1, sizeof (hio_error_stack_item_t));
  if (NULL == new_item) {
    /* not much can be done here. we are just plain OOM. */
    return;
  }

  new_item->hrc = hio_err_mpi(mpirc);

  /* TODO -- Should probably do somthing smarter here */
  new_item->error_string = malloc (strlen (temp) + 3 + resultlen);
  if (NULL == temp) {
    free (new_item);
    free (temp);
    return;
  }

  /* append the mpi error to the hio error string */
  strcpy (new_item->error_string, temp);
  strcat (new_item->error_string, ": ");
  strcat (new_item->error_string, mpi_error);

  /* done with this now */
  free (temp);

  /* push the error message onto the stack */
  if (NULL == context) {
    pthread_mutex_lock (&hio_error_stack_mutex);
    new_item->next = hio_error_stack_head;
    hio_error_stack_head = new_item;
    pthread_mutex_unlock (&hio_error_stack_mutex);
  } else {
    pthread_mutex_lock (&context->context_lock);
    new_item->next = (hio_error_stack_item_t *) context->context_error_stack;
    context->context_error_stack = (void *) new_item;
    pthread_mutex_unlock (&context->context_lock);
  }
}

int hio_err_mpi (int mpirc) {
  /* TODO: implement this */
  if (MPI_SUCCESS == mpirc) {
    return HIO_SUCCESS;
  }

  return HIO_ERROR;
}
#endif

int hio_err_get_last (hio_context_t context, char **error) {
  hio_error_stack_item_t *stack_error;
  int hrc;

  if (NULL == context) {
    pthread_mutex_lock (&hio_error_stack_mutex);
    stack_error = hio_error_stack_head;
    if (NULL != stack_error) {
      hio_error_stack_head = stack_error->next;
    }
    pthread_mutex_unlock (&hio_error_stack_mutex);
  } else {
    pthread_mutex_lock (&context->context_lock);
    stack_error = (hio_error_stack_item_t *) context->context_error_stack;
    if (NULL != stack_error) {
      context->context_error_stack = (void *) stack_error->next;
    }
    pthread_mutex_unlock (&context->context_lock);
  }

  if (NULL == stack_error) {
    /* no error */
    *error = NULL;
    return HIO_SUCCESS;
  }

  *error = stack_error->error_string;
  hrc = stack_error->hrc;
  free (stack_error);

  return hrc;
}

static int hio_err_print_last_vargs (hio_context_t context, FILE *output, char *format, va_list vargs) {
  char hostname[256] = "unknown";
  char datetime[30] = "unknown\n";
  char *hio_error;
  time_t timeval;
  int hrc, rc;

  /* dequeue the last error */
  hrc = hio_err_get_last (context, &hio_error);
  if (NULL == hio_error) {
    return 0;
  }

  /* try to get the hostname */
  (void) gethostname (hostname, 256);

  /* try to get the time */
  timeval = time (NULL);
  (void) ctime_r (&timeval, datetime);

  /* remove newline */
  datetime[strlen(datetime) - 1] = '\0';

  /* NTH: the following code prints a series of messages to the specified output
   * file handle. the code as is will probably not work properly if this function
   * is being called from multiple threads. in a future update this code should
   * be updated to buffer the error message before printing it out to the file
   * handle. */

  /* print out the timestamp */
  if (NULL == context) {
    rc = fprintf (output, "HIO %s <%s>: error code (%d) ", hostname, datetime, hrc);
  } else {
    rc = fprintf (output, "HIO %s <%s>: error code (%d) context (%s) ", hostname, datetime,
                  hrc, context->context_object.identifier);
  }

  /* print the user's error message */
  rc += vfprintf (output, format, vargs);

  /* finally, print out the hio error message */
  rc += fprintf (output, ": %s\n", hio_error);

  /* free the error message */
  free (hio_error);

  return rc;
}

int hio_err_print_last (hio_context_t ctx, FILE *output, char *format, ...) {
  va_list vargs;
  int rc;

  va_start (vargs, format);
  rc = hio_err_print_last_vargs (ctx, output, format, vargs);
  va_end (vargs);

  return rc;
}

int hio_err_print_all (hio_context_t ctx, FILE *output, char *format, ...)
{
  va_list vargs;
  int rc;

  /* loop until all error messages have been printed */
  do {
    va_start (vargs, format);
    rc = hio_err_print_last_vargs (ctx, output, format, vargs);
    va_end (vargs);

    if (0 == rc) {
      break;
    }
  } while (1);

  return HIO_SUCCESS;
}