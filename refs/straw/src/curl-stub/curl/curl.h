/* Stand-in for <curl/curl.h>: straw includes libcurl unconditionally for
 * reading remote files. This build reads local files only, so every curl
 * call fails; nothing links against libcurl. */
#ifndef OMICSWEAVE_CURL_STUB_H
#define OMICSWEAVE_CURL_STUB_H

typedef void CURL;
typedef int CURLcode;
typedef int CURLoption;

#define CURLE_OK 0
#define CURLOPT_URL 10002
#define CURLOPT_WRITEFUNCTION 20011
#define CURLOPT_WRITEDATA 10001
#define CURLOPT_RANGE 10007
#define CURLOPT_USERAGENT 10018
#define CURLOPT_FOLLOWLOCATION 52
#define CURLOPT_HEADERFUNCTION 20079

static inline CURL *curl_easy_init(void) { return 0; }
static inline CURLcode curl_easy_setopt(CURL *, CURLoption, ...) { return 1; }
static inline CURLcode curl_easy_perform(CURL *) { return 1; }
static inline void curl_easy_cleanup(CURL *) {}
static inline const char *curl_easy_strerror(CURLcode) { return "curl is not available in this build"; }

#endif
