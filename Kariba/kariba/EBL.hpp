#pragma once

#include <vector>

namespace kariba {

void ebl_atten_gil(size_t size, const std::vector<double> &en, std::vector<double> &lum,
                   double redshift);

}
