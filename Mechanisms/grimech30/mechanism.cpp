#include "mechanism.H"
const int rmap[NUM_REACTIONS] = {
  49,  51,  53,  55,  56,  58,  62,  69,  70,  71,  73,  75,  82,  84,  94,
  130, 139, 146, 157, 173, 240, 288, 303, 311, 317, 319, 11,  184, 236, 0,
  1,   32,  33,  34,  35,  36,  38,  39,  40,  41,  42,  165, 166, 186, 204,
  211, 226, 229, 268, 302, 2,   3,   4,   5,   6,   7,   8,   9,   10,  12,
  13,  14,  15,  16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,
  28,  29,  30,  31,  37,  43,  44,  45,  46,  47,  48,  50,  52,  54,  57,
  59,  60,  61,  63,  64,  65,  66,  67,  68,  72,  74,  76,  77,  78,  79,
  80,  81,  83,  85,  86,  87,  88,  89,  90,  91,  92,  93,  95,  96,  97,
  98,  99,  100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112,
  113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127,
  128, 129, 131, 132, 133, 134, 135, 136, 137, 138, 140, 141, 142, 143, 144,
  145, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 158, 159, 160, 161,
  162, 163, 164, 167, 168, 169, 170, 171, 172, 174, 175, 176, 177, 178, 179,
  180, 181, 182, 183, 185, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196,
  197, 198, 199, 200, 201, 202, 203, 205, 206, 207, 208, 209, 210, 212, 213,
  214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 227, 228, 230,
  231, 232, 233, 234, 235, 237, 238, 239, 241, 242, 243, 244, 245, 246, 247,
  248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262,
  263, 264, 265, 266, 267, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278,
  279, 280, 281, 282, 283, 284, 285, 286, 287, 289, 290, 291, 292, 293, 294,
  295, 296, 297, 298, 299, 300, 301, 304, 305, 306, 307, 308, 309, 310, 312,
  313, 314, 315, 316, 318, 320, 321, 322, 323, 324};

// Returns 0-based map of reaction order
void
GET_RMAP(int* _rmap)
{
  for (int j = 0; j < NUM_REACTIONS; ++j) {
    _rmap[j] = rmap[j];
  }
}

// Returns a count of gas species in a gas reaction, and their indices
// and stoichiometric coefficients. (Eq 50)
void
CKINU(const int i, int& nspec, int ki[], int nu[])
{
  const int ns[NUM_GAS_REACTIONS] = {
    2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 4, 2, 2, 2, 2, 3, 4, 4, 3, 4, 4, 4, 3,
    4, 3, 4, 3, 4, 3, 3, 4, 3, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 3, 3, 3, 4, 3, 4,
    3, 4, 4, 4, 4, 4, 4, 3, 4, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 3, 4, 4, 4, 5, 4, 3, 4, 3, 3, 4, 4, 4, 5, 4, 4, 3, 4, 4, 3,
    4, 4, 4, 4, 4, 4, 4, 2, 3, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 3, 4,
    4, 3, 4, 4, 4, 4, 3, 4, 4, 3, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 4,
    4, 3, 4, 4, 3, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 5, 4, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 5, 4, 4, 4, 4, 3, 4, 4, 3, 4, 4, 4, 4, 5, 5, 4, 5,
    5, 5, 3, 3, 5, 5, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 3, 4, 3, 4, 4, 4, 5, 3};
  const int kiv[NUM_GAS_REACTIONS * 5] = {
    2,  3,  0,  0,  0,  1,  2,  4,  0,  0,  0,  2,  1,  4,  0,  6,  2,  3,  4,
    0,  7,  2,  6,  4,  0,  9,  2,  14, 1,  0,  10, 2,  1,  16, 0,  11, 2,  14,
    0,  0,  11, 2,  1,  16, 0,  12, 2,  17, 1,  0,  13, 2,  12, 4,  0,  14, 2,
    15, 0,  0,  16, 2,  14, 4,  0,  16, 2,  15, 1,  0,  17, 2,  16, 4,  0,  18,
    2,  17, 4,  0,  19, 2,  17, 4,  0,  20, 2,  18, 4,  0,  20, 2,  19, 4,  0,
    21, 2,  9,  14, 0,  22, 2,  1,  27, 0,  22, 2,  21, 4,  0,  22, 2,  10, 14,
    0,  23, 2,  28, 1,  0,  24, 2,  12, 16, 0,  25, 2,  17, 12, 0,  26, 2,  25,
    4,  0,  27, 2,  14, 1,  0,  28, 2,  27, 4,  0,  28, 2,  10, 15, 0,  14, 3,
    15, 2,  0,  17, 3,  16, 6,  0,  1,  3,  6,  0,  0,  1,  3,  6,  0,  0,  1,
    3,  6,  0,  0,  1,  3,  6,  0,  0,  1,  3,  6,  0,  0,  1,  3,  2,  4,  0,
    1,  0,  0,  0,  0,  1,  0,  0,  0,  0,  1,  0,  0,  0,  0,  1,  0,  0,  0,
    0,  1,  4,  5,  0,  0,  1,  6,  5,  2,  0,  1,  6,  0,  3,  0,  1,  6,  4,
    0,  0,  1,  7,  0,  6,  0,  1,  7,  5,  4,  0,  9,  1,  8,  0,  0,  10, 1,
    12, 0,  0,  11, 1,  9,  0,  0,  12, 1,  13, 0,  0,  13, 1,  12, 0,  0,  1,
    16, 17, 0,  0,  1,  16, 14, 0,  0,  17, 1,  18, 0,  0,  17, 1,  19, 0,  0,
    17, 1,  0,  16, 0,  18, 1,  20, 0,  0,  18, 1,  17, 0,  0,  18, 1,  12, 4,
    0,  18, 1,  11, 5,  0,  19, 1,  20, 0,  0,  19, 1,  18, 1,  0,  19, 1,  17,
    0,  0,  19, 1,  12, 4,  0,  19, 1,  11, 5,  0,  20, 1,  18, 0,  0,  20, 1,
    19, 0,  0,  21, 1,  22, 0,  0,  22, 1,  23, 0,  0,  23, 1,  24, 0,  0,  23,
    1,  22, 0,  0,  24, 1,  25, 0,  0,  24, 1,  23, 0,  0,  25, 1,  26, 0,  0,
    25, 1,  24, 0,  0,  26, 1,  25, 0,  0,  1,  27, 11, 14, 0,  28, 1,  0,  27,
    0,  28, 1,  12, 14, 0,  1,  29, 28, 1,  0,  14, 0,  17, 0,  0,  0,  4,  1,
    5,  0,  4,  7,  0,  0,  0,  4,  5,  2,  0,  0,  6,  4,  5,  3,  0,  7,  4,
    5,  6,  0,  7,  4,  5,  6,  0,  8,  4,  14, 1,  0,  9,  4,  1,  16, 0,  10,
    4,  17, 1,  0,  10, 4,  9,  5,  0,  11, 4,  17, 1,  0,  12, 4,  20, 0,  0,
    12, 4,  10, 5,  0,  12, 4,  11, 5,  0,  13, 4,  12, 5,  0,  14, 4,  15, 1,
    0,  16, 4,  14, 5,  0,  17, 4,  5,  16, 0,  18, 4,  17, 5,  0,  19, 4,  17,
    5,  0,  20, 4,  18, 5,  0,  20, 4,  19, 5,  0,  21, 4,  1,  27, 0,  22, 4,
    28, 1,  0,  22, 4,  1,  29, 0,  22, 4,  21, 5,  0,  22, 4,  12, 14, 0,  23,
    4,  22, 5,  0,  24, 4,  23, 5,  0,  26, 4,  25, 5,  0,  28, 4,  5,  27, 0,
    6,  7,  3,  0,  0,  6,  7,  3,  0,  0,  10, 6,  17, 4,  0,  12, 6,  13, 3,
    0,  12, 6,  19, 4,  0,  14, 6,  15, 4,  0,  17, 6,  7,  16, 0,  8,  3,  14,
    2,  0,  8,  10, 21, 1,  0,  8,  12, 22, 1,  0,  9,  3,  16, 2,  0,  9,  0,
    10, 1,  0,  9,  5,  17, 1,  0,  9,  10, 22, 1,  0,  9,  12, 23, 1,  0,  9,
    13, 24, 1,  0,  9,  14, 27, 0,  0,  9,  15, 14, 16, 0,  9,  17, 28, 1,  0,
    9,  27, 22, 14, 0,  10, 3,  14, 1,  4,  10, 0,  12, 1,  0,  10, 22, 0,  0,
    0,  10, 12, 24, 1,  0,  10, 13, 12, 0,  0,  10, 14, 28, 0,  0,  10, 27, 23,
    14, 0,  11, 47, 10, 47, 0,  48, 11, 48, 10, 0,  11, 3,  14, 1,  4,  11, 3,
    14, 5,  0,  11, 0,  12, 1,  0,  11, 5,  20, 0,  0,  11, 5,  10, 5,  0,  11,
    12, 24, 1,  0,  11, 13, 12, 0,  0,  11, 14, 10, 14, 0,  11, 15, 10, 15, 0,
    11, 15, 17, 14, 0,  26, 11, 25, 12, 0,  12, 3,  19, 2,  0,  12, 3,  17, 4,
    0,  12, 7,  13, 6,  0,  12, 26, 0,  0,  0,  12, 25, 1,  0,  0,  12, 16, 13,
    14, 0,  17, 12, 13, 16, 0,  12, 20, 18, 13, 0,  12, 20, 19, 13, 0,  24, 12,
    23, 13, 0,  26, 12, 25, 13, 0,  16, 14, 1,  0,  0,  16, 14, 1,  0,  0,  16,
    3,  14, 6,  0,  18, 3,  17, 6,  0,  19, 3,  17, 6,  0,  21, 3,  14, 16, 0,
    21, 0,  22, 1,  0,  23, 3,  17, 16, 0,  24, 22, 0,  0,  0,  25, 3,  24, 6,
    0,  27, 3,  14, 4,  0,  27, 22, 14, 0,  0,  30, 35, 47, 2,  0,  30, 3,  35,
    2,  0,  30, 4,  1,  35, 0,  37, 2,  47, 3,  0,  37, 2,  35, 0,  0,  1,  37,
    47, 4,  0,  37, 4,  6,  47, 0,  37, 47, 2,  0,  0,  6,  35, 36, 4,  0,  35,
    2,  36, 0,  0,  36, 2,  35, 3,  0,  1,  36, 35, 4,  0,  31, 2,  1,  35, 0,
    1,  31, 0,  30, 0,  31, 4,  1,  38, 0,  31, 4,  5,  30, 0,  31, 3,  38, 2,
    0,  31, 3,  35, 4,  0,  30, 31, 1,  47, 0,  5,  31, 0,  38, 0,  31, 35, 47,
    4,  0,  31, 35, 1,  37, 0,  32, 2,  31, 4,  0,  32, 2,  1,  38, 0,  1,  32,
    0,  31, 0,  32, 4,  5,  31, 0,  34, 1,  47, 0,  0,  34, 1,  47, 0,  0,  34,
    3,  6,  47, 0,  34, 2,  47, 4,  0,  34, 2,  31, 35, 0,  1,  34, 0,  47, 0,
    34, 4,  5,  47, 0,  12, 34, 13, 47, 0,  1,  35, 38, 0,  0,  38, 2,  35, 4,
    0,  1,  38, 0,  35, 0,  38, 4,  5,  35, 0,  38, 3,  6,  35, 0,  39, 2,  14,
    30, 0,  39, 4,  1,  46, 0,  39, 5,  40, 4,  0,  39, 3,  46, 2,  0,  39, 0,
    1,  40, 0,  46, 2,  14, 35, 0,  1,  46, 14, 31, 0,  46, 4,  14, 1,  35, 30,
    46, 14, 47, 0,  46, 3,  15, 35, 0,  46, 14, 30, 0,  0,  46, 35, 14, 37, 0,
    46, 35, 15, 47, 0,  40, 39, 1,  0,  0,  40, 2,  1,  46, 0,  40, 2,  14, 31,
    0,  40, 2,  39, 4,  0,  40, 4,  1,  44, 0,  40, 4,  1,  45, 0,  40, 4,  14,
    32, 0,  1,  40, 41, 0,  0,  41, 30, 10, 47, 0,  8,  47, 39, 30, 0,  9,  47,
    40, 30, 0,  9,  47, 42, 0,  0,  10, 47, 40, 31, 0,  11, 47, 40, 31, 0,  8,
    35, 39, 2,  0,  8,  35, 14, 30, 0,  9,  35, 40, 2,  0,  9,  35, 1,  46, 0,
    9,  35, 16, 30, 0,  10, 35, 1,  45, 0,  10, 35, 40, 4,  0,  10, 35, 1,  43,
    0,  11, 35, 1,  45, 0,  11, 35, 40, 4,  0,  11, 35, 1,  43, 0,  12, 35, 5,
    40, 0,  12, 35, 41, 4,  0,  42, 2,  14, 1,  47, 42, 2,  40, 35, 0,  42, 3,
    16, 47, 2,  42, 4,  1,  16, 47, 1,  42, 10, 47, 0,  45, 2,  15, 31, 0,  45,
    2,  14, 38, 0,  45, 2,  46, 4,  0,  1,  45, 14, 32, 0,  1,  45, 0,  46, 0,
    45, 4,  5,  46, 0,  45, 4,  15, 32, 0,  45, 14, 31, 0,  0,  1,  43, 1,  45,
    0,  1,  43, 40, 4,  0,  1,  43, 14, 32, 0,  1,  44, 1,  45, 0,  27, 35, 14,
    43, 0,  12, 30, 1,  41, 0,  12, 30, 0,  40, 0,  1,  33, 0,  32, 0,  33, 4,
    5,  32, 0,  33, 2,  32, 4,  0,  15, 31, 14, 38, 0,  39, 36, 46, 35, 0,  46,
    36, 15, 37, 0,  15, 30, 14, 35, 0,  12, 2,  14, 1,  0,  24, 2,  51, 1,  0,
    25, 2,  52, 1,  0,  6,  4,  5,  3,  0,  12, 4,  17, 0,  0,  9,  0,  12, 0,
    0,  10, 3,  15, 1,  0,  10, 3,  17, 2,  0,  10, 22, 1,  0,  0,  11, 5,  17,
    0,  0,  23, 3,  51, 2,  0,  23, 3,  22, 6,  0,  52, 2,  51, 4,  0,  52, 2,
    12, 14, 4,  52, 3,  12, 14, 6,  52, 1,  51, 0,  0,  52, 1,  12, 14, 0,  52,
    4,  12, 14, 5,  52, 6,  12, 14, 7,  52, 13, 14, 0,  0,  28, 1,  51, 0,  0,
    51, 2,  10, 15, 1,  51, 3,  17, 14, 4,  51, 3,  16, 4,  0,  51, 1,  12, 16,
    0,  51, 1,  28, 0,  0,  51, 4,  28, 5,  0,  51, 4,  18, 16, 0,  25, 12, 50,
    0,  0,  50, 2,  49, 4,  0,  50, 1,  49, 0,  0,  50, 4,  49, 5,  0,  49, 7,
    50, 6,  0,  50, 12, 49, 13, 0,  24, 12, 49, 0,  0,  49, 2,  25, 17, 0,  49,
    1,  50, 0,  0,  49, 1,  25, 12, 0,  49, 4,  25, 18, 0,  49, 6,  50, 3,  0,
    49, 6,  25, 17, 4,  49, 12, 25, 0,  0};
  const int nuv[NUM_GAS_REACTIONS * 5] = {
    -2, 1,  0, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 2, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -2, 1,  0, 0, 0, -2, 1,  0, 0, 0,
    -2, 1,  0, 0, 0, -2, 1,  0, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 2, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0,
    -2, 1,  0, 0, 0, -2, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -2, 1,  1, 0, 0, -2, 1,  1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 1, -1, -1, 1, 1, 0,
    -2, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 2, 0, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 1,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 2, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -2, 1,  0, 0, 0, -2, 1,  1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, 1,  1, 0, 0, -1, 1,  1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 2, 1, 0,
    -2, 1,  2, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 2, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, 1,  1, 0, 0,
    -1, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 1,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, 1,  1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 1, -1, -1, 1, 1, 0, -1, -1, 1, 1, 1, -1, -1, 1, 1, 1,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 1,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 2, 0, -1, -1, 1, 1, 0, -2, 1,  2, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 1, -1, -1, 1, 1, 1, -1, -1, 1, 1, 0, -1, -1, 1, 1, 1,
    -1, -1, 1, 1, 1, -1, -1, 1, 1, 1, -1, 1,  1, 0, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 1, 1, -1, -1, 1, 1, 1, -1, -1, 2, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 1,
    -1, -1, 2, 0, 0};
  if (i < 1) {
    // Return max num species per reaction
    nspec = 5;
  } else {
    if (i > NUM_GAS_REACTIONS) {
      nspec = -1;
    } else {
      nspec = ns[i - 1];
      for (int j = 0; j < nspec; ++j) {
        ki[j] = kiv[(i - 1) * 5 + j] + 1;
        nu[j] = nuv[(i - 1) * 5 + j];
      }
    }
  }
}

// Returns the progress rates of each reactions
// Given P, T, and mole fractions
void
CKKFKR(
  const amrex::Real P,
  const amrex::Real T,
  const amrex::Real x[],
  amrex::Real q_f[],
  amrex::Real q_r[])
{
  amrex::Real c[53]; // temporary storage
  amrex::Real PORT =
    1e6 * P / (8.31446261815324e+07 * T); // 1e6 * P/RT so c goes to SI units

  // Compute conversion, see Eq 10
  for (int id = 0; id < 53; ++id) {
    c[id] = x[id] * PORT;
  }

  // convert to chemkin units
  progressRateFR(q_f, q_r, c, T);

  // convert to chemkin units
  for (int id = 0; id < 325; ++id) {
    q_f[id] *= 1.0e-6;
    q_r[id] *= 1.0e-6;
  }
}

// compute the progress rate for each reaction
// USES progressRate : todo switch to GPU
void
progressRateFR(
  amrex::Real* q_f, amrex::Real* q_r, amrex::Real* sc, amrex::Real T)
{
  const amrex::Real invT = 1.0 / T;
  const amrex::Real logT = log(T);
  // compute the Gibbs free energy
  amrex::Real g_RT[53];
  gibbs(g_RT, T);

  amrex::Real sc_qss[1];
  comp_qfqr(q_f, q_r, sc, sc_qss, T, invT, logT);
}

// save atomic weights into array
void
atomicWeight(amrex::Real* awt)
{
  awt[0] = 15.999000; // O
  awt[1] = 1.008000;  // H
  awt[2] = 12.011000; // C
  awt[3] = 14.007000; // N
  awt[4] = 39.950000; // Ar
}

// get atomic weight for all elements
void
CKAWT(amrex::Real* awt)
{
  atomicWeight(awt);
}

// Returns the elemental composition
// of the speciesi (mdim is num of elements)
void
CKNCF(int* ncf)
{
  int kd = 5;
  // Zero ncf
  for (int id = 0; id < kd * 53; ++id) {
    ncf[id] = 0;
  }

  // H2
  ncf[0 * kd + 1] = 2; // H

  // H
  ncf[1 * kd + 1] = 1; // H

  // O
  ncf[2 * kd + 0] = 1; // O

  // O2
  ncf[3 * kd + 0] = 2; // O

  // OH
  ncf[4 * kd + 1] = 1; // H
  ncf[4 * kd + 0] = 1; // O

  // H2O
  ncf[5 * kd + 1] = 2; // H
  ncf[5 * kd + 0] = 1; // O

  // HO2
  ncf[6 * kd + 1] = 1; // H
  ncf[6 * kd + 0] = 2; // O

  // H2O2
  ncf[7 * kd + 1] = 2; // H
  ncf[7 * kd + 0] = 2; // O

  // C
  ncf[8 * kd + 2] = 1; // C

  // CH
  ncf[9 * kd + 2] = 1; // C
  ncf[9 * kd + 1] = 1; // H

  // CH2
  ncf[10 * kd + 2] = 1; // C
  ncf[10 * kd + 1] = 2; // H

  // CH2(S)
  ncf[11 * kd + 2] = 1; // C
  ncf[11 * kd + 1] = 2; // H

  // CH3
  ncf[12 * kd + 2] = 1; // C
  ncf[12 * kd + 1] = 3; // H

  // CH4
  ncf[13 * kd + 2] = 1; // C
  ncf[13 * kd + 1] = 4; // H

  // CO
  ncf[14 * kd + 2] = 1; // C
  ncf[14 * kd + 0] = 1; // O

  // CO2
  ncf[15 * kd + 2] = 1; // C
  ncf[15 * kd + 0] = 2; // O

  // HCO
  ncf[16 * kd + 2] = 1; // C
  ncf[16 * kd + 1] = 1; // H
  ncf[16 * kd + 0] = 1; // O

  // CH2O
  ncf[17 * kd + 2] = 1; // C
  ncf[17 * kd + 1] = 2; // H
  ncf[17 * kd + 0] = 1; // O

  // CH2OH
  ncf[18 * kd + 2] = 1; // C
  ncf[18 * kd + 1] = 3; // H
  ncf[18 * kd + 0] = 1; // O

  // CH3O
  ncf[19 * kd + 2] = 1; // C
  ncf[19 * kd + 1] = 3; // H
  ncf[19 * kd + 0] = 1; // O

  // CH3OH
  ncf[20 * kd + 2] = 1; // C
  ncf[20 * kd + 1] = 4; // H
  ncf[20 * kd + 0] = 1; // O

  // C2H
  ncf[21 * kd + 2] = 2; // C
  ncf[21 * kd + 1] = 1; // H

  // C2H2
  ncf[22 * kd + 2] = 2; // C
  ncf[22 * kd + 1] = 2; // H

  // C2H3
  ncf[23 * kd + 2] = 2; // C
  ncf[23 * kd + 1] = 3; // H

  // C2H4
  ncf[24 * kd + 2] = 2; // C
  ncf[24 * kd + 1] = 4; // H

  // C2H5
  ncf[25 * kd + 2] = 2; // C
  ncf[25 * kd + 1] = 5; // H

  // C2H6
  ncf[26 * kd + 2] = 2; // C
  ncf[26 * kd + 1] = 6; // H

  // HCCO
  ncf[27 * kd + 2] = 2; // C
  ncf[27 * kd + 1] = 1; // H
  ncf[27 * kd + 0] = 1; // O

  // CH2CO
  ncf[28 * kd + 2] = 2; // C
  ncf[28 * kd + 1] = 2; // H
  ncf[28 * kd + 0] = 1; // O

  // HCCOH
  ncf[29 * kd + 2] = 2; // C
  ncf[29 * kd + 1] = 2; // H
  ncf[29 * kd + 0] = 1; // O

  // N
  ncf[30 * kd + 3] = 1; // N

  // NH
  ncf[31 * kd + 1] = 1; // H
  ncf[31 * kd + 3] = 1; // N

  // NH2
  ncf[32 * kd + 1] = 2; // H
  ncf[32 * kd + 3] = 1; // N

  // NH3
  ncf[33 * kd + 1] = 3; // H
  ncf[33 * kd + 3] = 1; // N

  // NNH
  ncf[34 * kd + 1] = 1; // H
  ncf[34 * kd + 3] = 2; // N

  // NO
  ncf[35 * kd + 3] = 1; // N
  ncf[35 * kd + 0] = 1; // O

  // NO2
  ncf[36 * kd + 3] = 1; // N
  ncf[36 * kd + 0] = 2; // O

  // N2O
  ncf[37 * kd + 3] = 2; // N
  ncf[37 * kd + 0] = 1; // O

  // HNO
  ncf[38 * kd + 1] = 1; // H
  ncf[38 * kd + 3] = 1; // N
  ncf[38 * kd + 0] = 1; // O

  // CN
  ncf[39 * kd + 2] = 1; // C
  ncf[39 * kd + 3] = 1; // N

  // HCN
  ncf[40 * kd + 2] = 1; // C
  ncf[40 * kd + 1] = 1; // H
  ncf[40 * kd + 3] = 1; // N

  // H2CN
  ncf[41 * kd + 2] = 1; // C
  ncf[41 * kd + 1] = 2; // H
  ncf[41 * kd + 3] = 1; // N

  // HCNN
  ncf[42 * kd + 2] = 1; // C
  ncf[42 * kd + 1] = 1; // H
  ncf[42 * kd + 3] = 2; // N

  // HCNO
  ncf[43 * kd + 2] = 1; // C
  ncf[43 * kd + 1] = 1; // H
  ncf[43 * kd + 3] = 1; // N
  ncf[43 * kd + 0] = 1; // O

  // HOCN
  ncf[44 * kd + 2] = 1; // C
  ncf[44 * kd + 1] = 1; // H
  ncf[44 * kd + 3] = 1; // N
  ncf[44 * kd + 0] = 1; // O

  // HNCO
  ncf[45 * kd + 2] = 1; // C
  ncf[45 * kd + 1] = 1; // H
  ncf[45 * kd + 3] = 1; // N
  ncf[45 * kd + 0] = 1; // O

  // NCO
  ncf[46 * kd + 2] = 1; // C
  ncf[46 * kd + 3] = 1; // N
  ncf[46 * kd + 0] = 1; // O

  // N2
  ncf[47 * kd + 3] = 2; // N

  // AR
  ncf[48 * kd + 4] = 1; // Ar

  // C3H7
  ncf[49 * kd + 2] = 3; // C
  ncf[49 * kd + 1] = 7; // H

  // C3H8
  ncf[50 * kd + 2] = 3; // C
  ncf[50 * kd + 1] = 8; // H

  // CH2CHO
  ncf[51 * kd + 2] = 2; // C
  ncf[51 * kd + 1] = 3; // H
  ncf[51 * kd + 0] = 1; // O

  // CH3CHO
  ncf[52 * kd + 2] = 2; // C
  ncf[52 * kd + 1] = 4; // H
  ncf[52 * kd + 0] = 1; // O
}

// Returns the vector of strings of element names
void
CKSYME_STR(amrex::Vector<std::string>& ename)
{
  ename.resize(5);
  ename[0] = "O";
  ename[1] = "H";
  ename[2] = "C";
  ename[3] = "N";
  ename[4] = "Ar";
}

// Returns the vector of strings of species names
void
CKSYMS_STR(amrex::Vector<std::string>& kname)
{
  kname.resize(53);
  kname[0] = "H2";
  kname[1] = "H";
  kname[2] = "O";
  kname[3] = "O2";
  kname[4] = "OH";
  kname[5] = "H2O";
  kname[6] = "HO2";
  kname[7] = "H2O2";
  kname[8] = "C";
  kname[9] = "CH";
  kname[10] = "CH2";
  kname[11] = "CH2(S)";
  kname[12] = "CH3";
  kname[13] = "CH4";
  kname[14] = "CO";
  kname[15] = "CO2";
  kname[16] = "HCO";
  kname[17] = "CH2O";
  kname[18] = "CH2OH";
  kname[19] = "CH3O";
  kname[20] = "CH3OH";
  kname[21] = "C2H";
  kname[22] = "C2H2";
  kname[23] = "C2H3";
  kname[24] = "C2H4";
  kname[25] = "C2H5";
  kname[26] = "C2H6";
  kname[27] = "HCCO";
  kname[28] = "CH2CO";
  kname[29] = "HCCOH";
  kname[30] = "N";
  kname[31] = "NH";
  kname[32] = "NH2";
  kname[33] = "NH3";
  kname[34] = "NNH";
  kname[35] = "NO";
  kname[36] = "NO2";
  kname[37] = "N2O";
  kname[38] = "HNO";
  kname[39] = "CN";
  kname[40] = "HCN";
  kname[41] = "H2CN";
  kname[42] = "HCNN";
  kname[43] = "HCNO";
  kname[44] = "HOCN";
  kname[45] = "HNCO";
  kname[46] = "NCO";
  kname[47] = "N2";
  kname[48] = "AR";
  kname[49] = "C3H7";
  kname[50] = "C3H8";
  kname[51] = "CH2CHO";
  kname[52] = "CH3CHO";
}

// compute the sparsity pattern of the chemistry Jacobian
void
SPARSITY_INFO(int* nJdata, const int* consP, int NCELLS)
{
  amrex::GpuArray<amrex::Real, 2916> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 53> conc = {0.0};
  for (int n = 0; n < 53; n++) {
    conc[n] = 1.0 / 53.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 54; k++) {
    for (int l = 0; l < 54; l++) {
      if (Jac[54 * k + l] != 0.0) {
        nJdata_tmp = nJdata_tmp + 1;
      }
    }
  }

  *nJdata = NCELLS * nJdata_tmp;
}

// compute the sparsity pattern of the system Jacobian
void
SPARSITY_INFO_SYST(int* nJdata, const int* consP, int NCELLS)
{
  amrex::GpuArray<amrex::Real, 2916> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 53> conc = {0.0};
  for (int n = 0; n < 53; n++) {
    conc[n] = 1.0 / 53.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 54; k++) {
    for (int l = 0; l < 54; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[54 * k + l] != 0.0) {
          nJdata_tmp = nJdata_tmp + 1;
        }
      }
    }
  }

  *nJdata = NCELLS * nJdata_tmp;
}

// compute the sparsity pattern of the simplified (for preconditioning) system
// Jacobian
void
SPARSITY_INFO_SYST_SIMPLIFIED(int* nJdata, const int* consP)
{
  amrex::GpuArray<amrex::Real, 2916> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 53> conc = {0.0};
  for (int n = 0; n < 53; n++) {
    conc[n] = 1.0 / 53.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 54; k++) {
    for (int l = 0; l < 54; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[54 * k + l] != 0.0) {
          nJdata_tmp = nJdata_tmp + 1;
        }
      }
    }
  }

  nJdata[0] = nJdata_tmp;
}

// compute the sparsity pattern of the chemistry Jacobian in CSC format -- base
// 0
void
SPARSITY_PREPROC_CSC(int* rowVals, int* colPtrs, const int* consP, int NCELLS)
{
  amrex::GpuArray<amrex::Real, 2916> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 53> conc = {0.0};
  for (int n = 0; n < 53; n++) {
    conc[n] = 1.0 / 53.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int nc = 0; nc < NCELLS; nc++) {
    int offset_row = nc * 54;
    int offset_col = nc * 54;
    for (int k = 0; k < 54; k++) {
      for (int l = 0; l < 54; l++) {
        if (Jac[54 * k + l] != 0.0) {
          rowVals[nJdata_tmp] = l + offset_row;
          nJdata_tmp = nJdata_tmp + 1;
        }
      }
      colPtrs[offset_col + (k + 1)] = nJdata_tmp;
    }
  }
}

// compute the sparsity pattern of the chemistry Jacobian in CSR format -- base
// 0
void
SPARSITY_PREPROC_CSR(
  int* colVals, int* rowPtrs, const int* consP, int NCELLS, int base)
{
  amrex::GpuArray<amrex::Real, 2916> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 53> conc = {0.0};
  for (int n = 0; n < 53; n++) {
    conc[n] = 1.0 / 53.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtrs[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 54;
      for (int l = 0; l < 54; l++) {
        for (int k = 0; k < 54; k++) {
          if (Jac[54 * k + l] != 0.0) {
            colVals[nJdata_tmp - 1] = k + 1 + offset;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
        rowPtrs[offset + (l + 1)] = nJdata_tmp;
      }
    }
  } else {
    rowPtrs[0] = 0;
    int nJdata_tmp = 0;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 54;
      for (int l = 0; l < 54; l++) {
        for (int k = 0; k < 54; k++) {
          if (Jac[54 * k + l] != 0.0) {
            colVals[nJdata_tmp] = k + offset;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
        rowPtrs[offset + (l + 1)] = nJdata_tmp;
      }
    }
  }
}

// compute the sparsity pattern of the system Jacobian
// CSR format BASE is user choice
void
SPARSITY_PREPROC_SYST_CSR(
  int* colVals, int* rowPtr, const int* consP, int NCELLS, int base)
{
  amrex::GpuArray<amrex::Real, 2916> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 53> conc = {0.0};
  for (int n = 0; n < 53; n++) {
    conc[n] = 1.0 / 53.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 54;
      for (int l = 0; l < 54; l++) {
        for (int k = 0; k < 54; k++) {
          if (k == l) {
            colVals[nJdata_tmp - 1] = l + 1 + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[54 * k + l] != 0.0) {
              colVals[nJdata_tmp - 1] = k + 1 + offset;
              nJdata_tmp = nJdata_tmp + 1;
            }
          }
        }
        rowPtr[offset + (l + 1)] = nJdata_tmp;
      }
    }
  } else {
    rowPtr[0] = 0;
    int nJdata_tmp = 0;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 54;
      for (int l = 0; l < 54; l++) {
        for (int k = 0; k < 54; k++) {
          if (k == l) {
            colVals[nJdata_tmp] = l + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[54 * k + l] != 0.0) {
              colVals[nJdata_tmp] = k + offset;
              nJdata_tmp = nJdata_tmp + 1;
            }
          }
        }
        rowPtr[offset + (l + 1)] = nJdata_tmp;
      }
    }
  }
}

// compute the sparsity pattern of the simplified (for precond) system Jacobian
// on CPU BASE 0
void
SPARSITY_PREPROC_SYST_SIMPLIFIED_CSC(
  int* rowVals, int* colPtrs, int* indx, const int* consP)
{
  amrex::GpuArray<amrex::Real, 2916> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 53> conc = {0.0};
  for (int n = 0; n < 53; n++) {
    conc[n] = 1.0 / 53.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int k = 0; k < 54; k++) {
    for (int l = 0; l < 54; l++) {
      if (k == l) {
        rowVals[nJdata_tmp] = l;
        indx[nJdata_tmp] = 54 * k + l;
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[54 * k + l] != 0.0) {
          rowVals[nJdata_tmp] = l;
          indx[nJdata_tmp] = 54 * k + l;
          nJdata_tmp = nJdata_tmp + 1;
        }
      }
    }
    colPtrs[k + 1] = nJdata_tmp;
  }
}

// compute the sparsity pattern of the simplified (for precond) system Jacobian
// CSR format BASE is under choice
void
SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(
  int* colVals, int* rowPtr, const int* consP, int base)
{
  amrex::GpuArray<amrex::Real, 2916> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 53> conc = {0.0};
  for (int n = 0; n < 53; n++) {
    conc[n] = 1.0 / 53.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int l = 0; l < 54; l++) {
      for (int k = 0; k < 54; k++) {
        if (k == l) {
          colVals[nJdata_tmp - 1] = l + 1;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[54 * k + l] != 0.0) {
            colVals[nJdata_tmp - 1] = k + 1;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
      }
      rowPtr[l + 1] = nJdata_tmp;
    }
  } else {
    rowPtr[0] = 0;
    int nJdata_tmp = 0;
    for (int l = 0; l < 54; l++) {
      for (int k = 0; k < 54; k++) {
        if (k == l) {
          colVals[nJdata_tmp] = l;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[54 * k + l] != 0.0) {
            colVals[nJdata_tmp] = k;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
      }
      rowPtr[l + 1] = nJdata_tmp;
    }
  }
}
