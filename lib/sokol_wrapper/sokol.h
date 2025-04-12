#if defined(__APPLE__)
#define SOKOL_METAL
#elif defined(__EMSCRIPTEN__)
#define SOKOL_GLES3
#elif defined(_WIN32)
#define SOKOL_WIN32_FORCE_MAIN
#define SOKOL_D3D11
#elif defined(__ANDROID__)
#define SOKOL_GLES3
#elif defined(__linux__) || defined(__unix__)
#define SOKOL_GLCORE
#else
#error "Unknown platform"
#endif

#if defined(SOKOL_D3D11)
#define ORIGIN_TOP_LEFT true
#else
#define ORIGIN_TOP_LEFT false
#endif

#define SOKOL_DEBUG
#define SOKOL_TRACE_HOOKS
#include <sokol_log.h>
#include <sokol_app.h>
#include <sokol_gfx.h>
#include <sokol_glue.h>
#include <sokol_time.h>