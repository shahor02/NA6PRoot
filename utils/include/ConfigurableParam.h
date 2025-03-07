// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

// first version 8/2018, Sandro Wenzel

#ifndef COMMON_SIMCONFIG_INCLUDE_SIMCONFIG_CONFIGURABLEPARAM_H_
#define COMMON_SIMCONFIG_INCLUDE_SIMCONFIG_CONFIGURABLEPARAM_H_

#include <vector>
#include <cassert>
#include <map>
#include <unordered_map>
#include <boost/property_tree/ptree_fwd.hpp>
#include <typeinfo>
#include <iostream>
#include <array>

class TFile;
class TRootIOCtor;
class TDataMember;

namespace na6p
{
namespace conf
{
// Base class for a configurable parameter.
//
// A configurable parameter (ConfigurableParameter) is a simple class, defining
// a few (pod) properties/members which are registered
// in a global (boost) property tree / structure.
//
// The features that we provide here are:
// *) Automatic translation from C++ data description to INI/JSON/XML
//    format via ROOT introspection and boost property trees and
//    the possibility to readably save the configuration
// *) Serialization/Deserialization into ROOT binary blobs (for the purpose
//    of writing/retrieving parameters from CCDB and to pass parameters along the processing chain)
// *) Automatic integration of sub-classes into a common configuration
// *) Be able to query properties from high level interfaces (just knowing
// *) Be able to set properties from high-level interfaces (and modifying the underlying
//    C++ object)
// *) Automatic ability to modify parameters from the command-line
// *) Keeping track of the provenance of individual parameter values; The user is able
//    to see whether is parameter is defaulted-code-value/coming-from-CCDB/coming-from-comandline
//
// Note that concrete parameter sub-classes **must** be implemented
// by inheriting from ConfigurableParamHelper and not from this class.
//
// ---------------------
// Example: To define a parameter class TPCGasParameters, one does the following:
//
// class TPCGasParamer : public ConfigurableParamHelper<TPCGasParameter>
// {
//   public:
//     double getGasDensity() const { return mGasDensity; }
//   private: // define properties AND default values
//     double mGasDensity = 1.23;
//     int mGasMaterialID = 1;
//
//     NA6PParamDef(TPCGasParameter, TPCGas); // a macro implementing some magic
// }
//
//
// We can now query the parameters in various ways
// - All parameter classes are singletons and we can say: TPCGasParameter::Instance().getGasDensity();
// - We can query by key (using classname + parameter name) from the global registry:
// -    ConfigurableParameter::getValueAs<double>("TPCGas", "mGasDensity");
//
// We can modify the parameters via the global registry together with an automatic syncing
// of the underlying C++ object:
// - ConfigurableParameter::setValue("TPCGas.mGasDensity", "0.5");
//
// - TPCGasParameter::Instance().getGasParameter() will now return 0.5;
//
// This feature allows to easily modify parameters at run-time via a textual representation
// (for example by giving strings on the command line)
//
// The collection of all parameter keys and values can be stored to a human/machine readable
// file
//  - ConfigurableParameter::writeJSON("thisconfiguration.json")

struct EnumLegalValues {
  std::vector<std::pair<std::string, int>> vvalues;

  bool isLegal(const std::string& value) const
  {
    for (auto& v : vvalues) {
      if (v.first == value) {
        return true;
      }
    }
    return false;
  }

  bool isLegal(int value) const
  {
    for (auto& v : vvalues) {
      if (v.second == value) {
        return true;
      }
    }
    return false;
  }

  std::string toString() const;
  int getIntValue(const std::string& value) const;
};

class EnumRegistry
{
 public:
  void add(const std::string& key, const TDataMember* dm);

  bool contains(const std::string& key) const
  {
    return entries.count(key) > 0;
  }

  std::string toString() const;

  const EnumLegalValues* operator[](const std::string& key) const
  {
    auto iter = entries.find(key);
    return iter != entries.end() ? &iter->second : nullptr;
  }

 private:
  std::unordered_map<std::string, EnumLegalValues> entries;
};

class ConfigurableParam
{
 public:
  enum EParamProvenance {
    kCODE /* from default code initialization */,
    kCCDB /* overwritten from CCDB */,
    kRT /* changed during runtime via API call setValue (for example command line) */,
    kRTF /* changed during runtime via updateFromFile */,
    /* can add more modes here */
  };

  enum class EParamUpdateStatus {
    Changed,   // param was successfully changed
    Unchanged, // param was not changed: new value is the same as previous
    Failed     // failed to update param
  };

  static std::string toString(EParamProvenance p)
  {
    static std::array<std::string, 4> names = {"CODE", "CCDB", "RT", "RTF"};
    return names[(int)p];
  }

  // get the name of the configurable Parameter
  virtual std::string getName() const = 0;

  // print the current keys and values to screen (optionally with provenance information)
  virtual void printKeyValues(bool showprov = true, bool useLogger = false) const = 0;

  // get a single size_t hash_value of this parameter (can be used as a checksum to see
  // if object changed or different)
  virtual size_t getHash() const = 0;

  // return the provenance of the member key
  virtual EParamProvenance getMemberProvenance(const std::string& key) const = 0;

  static EParamProvenance getProvenance(const std::string& key);

  static void printAllRegisteredParamNames();
  static void printAllKeyValuePairs(bool useLogger = false);

  static const std::string& getOutputDir() { return sOutputDir; }

  static void setOutputDir(const std::string& d) { sOutputDir = d; }

  static bool configFileExists(std::string const& filepath);

  // writes a human readable JSON file of all parameters
  static void writeJSON(std::string const& filename, std::string const& keyOnly = "");
  // writes a human readable INI file of all parameters
  static void writeINI(std::string const& filename, std::string const& keyOnly = "");

  // can be used instead of using API on concrete child classes
  template <typename T>
  static T getValueAs(std::string key)
  {
    return [](auto* tree, const std::string& key) -> T {
      if (!sIsFullyInitialized) {
        initialize();
      }
      return tree->template get<T>(key);
    }(sPtree, key);
  }

  template <typename T>
  static void setValue(std::string const& mainkey, std::string const& subkey, T x, EParamProvenance prov = kRT)
  {
    if (!sIsFullyInitialized) {
      initialize();
    }
    return [&subkey, &x, &mainkey, &prov](auto* tree) -> void {
      assert(tree);
      try {
        auto key = mainkey + "." + subkey;
        if (tree->template get_optional<std::string>(key).is_initialized()) {
          tree->put(key, x);
          auto changed = updateThroughStorageMap(mainkey, subkey, typeid(T), (void*)&x);
          if (changed != EParamUpdateStatus::Failed) {
            sValueProvenanceMap->find(key)->second = prov; // set to runtime
          }
        }
      } catch (std::exception const& e) {
        std::cerr << "Error in setValue (T) " << e.what() << "\n";
      }
    }(sPtree);
  }

  static void setProvenance(std::string const& mainkey, std::string const& subkey, EParamProvenance p)
  {
    if (!sIsFullyInitialized) {
      std::cerr << "setProvenance was called on non-initialized ConfigurableParam\n";
      return;
    }
    try {
      auto key = mainkey + "." + subkey;
      auto keyProv = sValueProvenanceMap->find(key);
      if (keyProv != sValueProvenanceMap->end()) {
        keyProv->second = p;
      }
    } catch (std::exception const& e) {
      std::cerr << "Error in setProvenance (T) " << e.what() << "\n";
    }
  }

  // specialized for std::string
  // which means that the type will be converted internally
  static void setValue(std::string const& key, std::string const& valuestring, EParamProvenance p = kRT);
  static void setEnumValue(const std::string&, const std::string&, EParamProvenance p = kRT);
  static void setArrayValue(const std::string&, const std::string&, EParamProvenance p = kRT);

  // update the storagemap from a vector of key/value pairs, calling setValue for each pair
  static void setValues(std::vector<std::pair<std::string, std::string>> const& keyValues, EParamProvenance p = kRT);

  // initializes the parameter database
  static void initialize();

  // create CCDB snapsnot
  static void toCCDB(std::string filename);
  // load from (CCDB) snapshot
  static void fromCCDB(std::string filename);

  // allows to provide a string of key-values from which to update
  // (certain) key-values
  // propagates changes down to each registered configuration
  // might be useful to get stuff from the command line
  static void updateFromString(std::string const&);

  // provide a path to a configuration file with ConfigurableParam key/values
  // If nonempty comma-separated paramsList is provided, only those params will
  // be updated, absence of data for any of requested params will lead to fatal
  static void updateFromFile(std::string const&, std::string const& paramsList = "", bool unchangedOnly = false);

  // interface for use from the CCDB API; allows to sync objects read from CCDB with the information
  // stored in the registry; modifies given object as well as registry
  virtual void syncCCDBandRegistry(void* obj) = 0;

 protected:
  // constructor is doing nothing else but
  // registering the concrete parameters
  ConfigurableParam();

  friend std::ostream& operator<<(std::ostream& out, const ConfigurableParam& me);

  static void initPropertyTree();
  static EParamUpdateStatus updateThroughStorageMap(std::string, std::string, std::type_info const&, void*);
  static EParamUpdateStatus updateThroughStorageMapWithConversion(std::string const&, std::string const&);

  virtual ~ConfigurableParam() = default;

  // fill property tree with the key-values from the sub-classes
  virtual void putKeyValues(boost::property_tree::ptree*) = 0;
  virtual void output(std::ostream& out) const = 0;

  virtual void serializeTo(TFile*) const = 0;
  virtual void initFrom(TFile*) = 0;

  // static map keeping, for each configuration key, its memory location and type
  // (internal use to easily sync updates, this is ok since parameter classes are singletons)
  static std::map<std::string, std::pair<std::type_info const&, void*>>* sKeyToStorageMap;

  // keep track of provenance of parameters and values
  static std::map<std::string, ConfigurableParam::EParamProvenance>* sValueProvenanceMap;

  // A registry of enum names and their allowed values
  // (stored as a vector of pairs <enumValueLabel, enumValueInt>)
  static EnumRegistry* sEnumRegistry;

  static std::string sOutputDir;

  void setRegisterMode(bool b) { sRegisterMode = b; }
  bool isInitialized() const { return sIsFullyInitialized; }

  // friend class o2::ccdb::CcdbApi;
 private:
  // static registry for implementations of this type
  static std::vector<ConfigurableParam*>* sRegisteredParamClasses; //!
  // static property tree (stocking all key - value pairs from instances of type ConfigurableParam)
  static boost::property_tree::ptree* sPtree; //!
  static bool sIsFullyInitialized;            //!
  static bool sRegisterMode;                  //! (flag to enable/disable autoregistering of child classes)
};

} // end namespace conf
} // namespace na6p

// a helper macro for boilerplate code in parameter classes
#define NA6PParamDef(classname, key)                \
 public:                                            \
  classname(TRootIOCtor*) {}                        \
  classname(classname const&) = delete;             \
                                                    \
 private:                                           \
  static constexpr char const* const sKey = key;    \
  static classname sInstance;                       \
  classname() = default;                            \
  template <typename T>                             \
  friend class na6p::conf::ConfigurableParamHelper; \
  template <typename T, typename P>                 \
  friend class na6p::conf::ConfigurableParamPromoter;

// a helper macro to implement necessary symbols in source
#define O2ParamImpl(classname) classname classname::sInstance;

#endif /* COMMON_SIMCONFIG_INCLUDE_SIMCONFIG_CONFIGURABLEPARAM_H_ */
