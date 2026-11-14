// C ABI over straw's C++ entry point: converts C strings in and the record
// vector out, nothing else. Results are released with straw_free.

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "straw.h"

extern "C" {

struct straw_result {
    int64_t n;
    int32_t *bin_x;
    int32_t *bin_y;
    float *counts;
};

// straw(matrix_type, norm, file, loc1, loc2, unit, bin_size); NULL on an
// exception.
straw_result *straw_query(const char *matrix_type, const char *norm, const char *file, const char *loc1,
                          const char *loc2, const char *unit, int32_t bin_size) {
    try {
        std::vector<contactRecord> records = straw(matrix_type, norm, file, loc1, loc2, unit, bin_size);
        auto *out = static_cast<straw_result *>(std::calloc(1, sizeof(straw_result)));
        out->n = static_cast<int64_t>(records.size());
        size_t n = records.size() ? records.size() : 1;
        out->bin_x = static_cast<int32_t *>(std::malloc(n * sizeof(int32_t)));
        out->bin_y = static_cast<int32_t *>(std::malloc(n * sizeof(int32_t)));
        out->counts = static_cast<float *>(std::malloc(n * sizeof(float)));
        for (size_t i = 0; i < records.size(); i++) {
            out->bin_x[i] = records[i].binX;
            out->bin_y[i] = records[i].binY;
            out->counts[i] = records[i].counts;
        }
        return out;
    } catch (...) {
        return nullptr;
    }
}

void straw_free(straw_result *r) {
    if (!r) return;
    std::free(r->bin_x);
    std::free(r->bin_y);
    std::free(r->counts);
    std::free(r);
}
}
