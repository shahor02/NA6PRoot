// NA6PCCopyright

#include <vector>
#include <array>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include "StringUtils.h"
#include "fairlogger/Logger.h"
#include "MagneticFieldRegion.h"

void MagneticFieldRegion::loadFlukaField(const std::string& filename)
{
  // read field provided by MGNCREAT / MGNDATA fluka cards
  if (mName.empty()) {
    mName = filename;
  }
  std::ifstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("Cannot open file: " + filename);
  }
  std::string line;
  int seenMGNCREAT = 0;
  std::vector<float> MGNCREATcache;
  while (std::getline(file, line)) {
    if (line.find("MGNCREAT") != std::string::npos) {
      seenMGNCREAT = true;
      cacheValues(line, MGNCREATcache);
    } else if (line.find("MGNDATA") != std::string::npos) {
      if (seenMGNCREAT == 1) {
        if (MGNCREATcache.size() >= 13) {
          mXMin = MGNCREATcache[7];
          mYMin = MGNCREATcache[8];
          mZMin = MGNCREATcache[9];
          mXMax = MGNCREATcache[10];
          mYMax = MGNCREATcache[11];
          mZMax = MGNCREATcache[12];
          mNX = static_cast<int>(MGNCREATcache[4]);
          mNY = static_cast<int>(MGNCREATcache[5]);
          mNZ = static_cast<int>(MGNCREATcache[6]);
          mNX1 = mNX - 1;
          mNY1 = mNY - 1;
          mNZ1 = mNZ - 1;
          mDX = (mXMax - mXMin) / mNX1;
          mDY = (mYMax - mYMin) / mNY1;
          mDZ = (mZMax - mZMin) / mNZ1;
          mDXI = 1. / mDX;
          mDYI = 1. / mDY;
          mDZI = 1. / mDZ;
          seenMGNCREAT = 2;
          mFieldData.reserve(3 * mNX * mNY * mNZ);
        } else {
          LOGP(fatal, "MGNCREAT collected only {} values", MGNCREATcache.size());
        }
      }
      cacheValues(line, mFieldData);
    }
  }
  file.close();
  if (mFieldData.size() != size_t(3 * mNX * mNY * mNZ)) {
    LOGP(fatal, "Expected {} field points for {}x{}x{} grid, loaded {}", 3 * mNX * mNY * mNZ, mNX, mNY, mNZ, mFieldData.size());
  }
}

void MagneticFieldRegion::loadOpera3DField(const std::string& filename)
{
  if (mName.empty()) {
    mName = filename;
  }

  std::ifstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("Cannot open file: " + filename);
  }

  std::string line;
  double xmin_orig, xmax_orig, ymin_orig, ymax_orig, zmin_orig, zmax_orig;
  int nx_orig, ny_orig, nz_orig;

  // Read header parameters (original grid)
  while (std::getline(file, line)) {
    if (line.empty() || line[0] == '!') {
      continue;
    }

    size_t pos = line.find('>');
    if (pos != std::string::npos) {
      std::string param = line.substr(0, pos);
      std::string value = line.substr(pos + 1);

      param.erase(0, param.find_first_not_of(" \t"));
      param.erase(param.find_last_not_of(" \t") + 1);
      value.erase(0, value.find_first_not_of(" \t"));
      value.erase(value.find_last_not_of(" \t") + 1);

      if (param == "xmin")
        xmin_orig = std::stod(value);
      else if (param == "xmax")
        xmax_orig = std::stod(value);
      else if (param == "ymin")
        ymin_orig = std::stod(value);
      else if (param == "ymax")
        ymax_orig = std::stod(value);
      else if (param == "zmin")
        zmin_orig = std::stod(value);
      else if (param == "zmax")
        zmax_orig = std::stod(value);
      else if (param == "nx")
        nx_orig = std::stoi(value);
      else if (param == "ny")
        ny_orig = std::stoi(value);
      else if (param == "nz")
        nz_orig = std::stoi(value);
    } else {
      // First data line - rewind
      std::streampos currentPos = file.tellg();
      file.seekg(currentPos - static_cast<std::streamoff>(line.length() + 1));
      break;
    }
  }

  // Read original field data into temporary storage
  std::vector<std::array<float, 3>> originalData;
  originalData.reserve(nx_orig * ny_orig * nz_orig);

  while (std::getline(file, line)) {
    if (line.empty() || line[0] == '!')
      continue;

    std::istringstream iss(line);
    double x, y, z, fx, fy, fz;

    if (iss >> x >> y >> z >> fx >> fy >> fz) {
      originalData.push_back({static_cast<float>(fx), static_cast<float>(fy), static_cast<float>(fz)});
    }
  }
  file.close();

  // Set up extended grid parameters (symmetric around origin)
  mXMin = -xmax_orig;
  mXMax = xmax_orig;
  mYMin = -ymax_orig;
  mYMax = ymax_orig;
  mZMin = -zmax_orig;
  mZMax = zmax_orig;

  // Extended grid dimensions (2n-1 to avoid duplicating center points)
  mNX = 2 * nx_orig - 1;
  mNY = 2 * ny_orig - 1;
  mNZ = 2 * nz_orig - 1;

  mNX1 = mNX - 1;
  mNY1 = mNY - 1;
  mNZ1 = mNZ - 1;

  mDX = (mXMax - mXMin) / mNX1;
  mDY = (mYMax - mYMin) / mNY1;
  mDZ = (mZMax - mZMin) / mNZ1;
  mDXI = 1. / mDX;
  mDYI = 1. / mDY;
  mDZI = 1. / mDZ;

  // Allocate extended field data
  mFieldData.clear();
  mFieldData.resize(3 * mNX * mNY * mNZ);

  // Helper lambda to set field data in extended grid
  auto setFieldData = [this](int ix, int iy, int iz, float Bx, float By, float Bz) {
    int index = (iz * mNY + iy) * mNX + ix;
    mFieldData[3 * index + 0] = Bx;
    mFieldData[3 * index + 1] = By;
    mFieldData[3 * index + 2] = Bz;
  };

  // Fill extended grid with symmetric data
  for (int ix = 0; ix < mNX; ++ix) {
    for (int iy = 0; iy < mNY; ++iy) {
      for (int iz = 0; iz < mNZ; ++iz) {
        // Map extended grid indices to original grid
        auto reflect_index = [](int i, int N) {
            return (i < N) ? (N - 1 - i) : (i - N + 1);
        };

        int ix_orig = reflect_index(ix, nx_orig);
        int iy_orig = reflect_index(iy, ny_orig);
        int iz_orig = reflect_index(iz, nz_orig);

        int origIndex = (iz_orig * ny_orig + iy_orig) * nx_orig + ix_orig;
        float Bx = originalData[origIndex][0];
        float By = originalData[origIndex][1];
        float Bz = originalData[origIndex][2];

        // Apply dipole symmetry transformations
        // Determine signs based on which quadrant we're in
        bool neg_x = (ix < nx_orig);
        bool neg_y = (iy < ny_orig);
        bool neg_z = (iz < nz_orig);

        // Dipole symmetry:

        // By(x,y,z) = By(-x,y,z) = By(x,-y,z) = By(x,y,-z)
        // Bx(x,y,z) = -Bx(-x,y,z) = Bx(x,-y,z) = Bx(x,y,-z)
        // Bz(x,y,z) = Bz(-x,y,z) = -Bz(x,-y,z) = -Bz(x,y,-z)

        if (neg_x) {
          Bx *= -1;
        } // Only Bx changes sign
        if (neg_y) {
          Bz *= -1;
        } // Only Bz changes sign
        if (neg_z) {
          Bz *= -1;
        } // Only Bz changes sign

        setFieldData(ix, iy, iz, Bx, By, Bz);
      }
    }
  }

  LOGP(info, "Loaded symmetric field from {} with extended grid {}x{}x{} (from original {}x{}x{})",
       filename, mNX, mNY, mNZ, nx_orig, ny_orig, nz_orig);
  LOGP(info, "Extended bounds: {}< X <{}, {}< Y <{}, {}< Z <{}",
       mXMin, mXMax, mYMin, mYMax, mZMin, mZMax);
}

void MagneticFieldRegion::cacheValues(const std::string& line, std::vector<float>& cachev)
{
  std::istringstream iss(line);
  auto vstr = na6p::utils::Str::tokenize(line, ',', true, false);
  for (auto& s : vstr) {
    try {
      cachev.push_back(std::stof(s));
    } catch (...) {
      // ignore
    }
  }
}

void MagneticFieldRegion::setRefPosition(float x, float y, float z)
{
  std::string rep;
  mRefPos[0] = x;
  mRefPos[1] = y;
  mRefPos[2] = z;
  mBoxPos[0][0] = mXMin + mRefPos[0];
  mBoxPos[0][1] = mXMax + mRefPos[0];
  rep += fmt::format("{:.1f}<X<{:.1f} ", mBoxPos[0][0], mBoxPos[0][1]);
  mBoxPos[1][0] = mYMin + mRefPos[1];
  mBoxPos[1][1] = mYMax + mRefPos[1];
  rep += fmt::format("{:.1f}<Y<{:.1f} ", mBoxPos[1][0], mBoxPos[1][1]);
  mBoxPos[2][0] = mZMin + mRefPos[2];
  mBoxPos[2][1] = mZMax + mRefPos[2];
  rep += fmt::format("{:.1f}<Z<{:.1f} ", mBoxPos[2][0], mBoxPos[2][1]);
  LOGP(info, "Adding Field {} in {}", mName, rep);
}
