#pragma once

#include <cstdlib>
#include <fstream>
#include <map>
#include <set>
#include <stdexcept>

// Enum used to label types of analysis ntuple files
enum class NtupleFileType {

  // Data taken with the beam on
  kOnBNB,

  // Data taken with the beam off
  kExtBNB,

  // MC events (standard BNB)
  kNumuMC,

  // MC events (intrinsic-nue-enhanced BNB)
  kIntrinsicNueMC,

  // MC events (dirt BNB)
  kDirtMC,

  // *** DetVar MC ***
  kDetVarMCCV, // central value
  kDetVarMCLYatten, // light-yield attenuation
  kDetVarMCLYdown, // light-yield down
  kDetVarMCLYrayl, // light-yield Rayleigh scattering
  kDetVarMCRecomb2, // light-yield recombination 2
  kDetVarMCSCE, // space charge effect
  kDetVarMCWMAngleXZ, // wireMod angle XZ
  kDetVarMCWMAngleYZ, // wireMod angle YZ
  kDetVarMCWMdEdx, // wireMod dE/dx
  kDetVarMCWMX, // wireMod X
  kDetVarMCWMYZ, // wireMod YZ
};

// Singleton class that keeps track of the various ntuple files to be analyzed
class FilePropertiesManager {

  public:

    // This is a singleton class, so we'll delete the copy constructor
    // the move constructor, and the assignment operators
    FilePropertiesManager( const FilePropertiesManager& ) = delete;
    FilePropertiesManager( FilePropertiesManager&& ) = delete;
    FilePropertiesManager& operator=( const FilePropertiesManager& )
      = delete;
    FilePropertiesManager& operator=( FilePropertiesManager&& )
      = delete;

    // Get a const reference to the singleton instance of the
    // FilePropertiesManager
    inline static const FilePropertiesManager& Instance() {

      // Create the FilePropertiesManager object using a static variable.
      // This ensures that the singleton instance is only created once.
      static std::unique_ptr<FilePropertiesManager>
        the_instance( new FilePropertiesManager() );

      // Return a reference to the singleton instance
      return *the_instance;
    }

    // Simple container that stores the number of triggers and the POT exposure
    // represented by a particular data ntuple
    struct TriggersAndPOT {
      TriggersAndPOT() {}
      TriggersAndPOT( int trig, double pot )
        : trigger_count_( trig ), pot_( pot ) {}

      int trigger_count_ = 0;
      double pot_ = 0.;
    };

    inline const std::map< std::string, TriggersAndPOT >& data_norm_map() const
      { return data_norm_map_; }


    inline const std::map< int, std::map<NtupleFileType,
      std::set<std::string> > >& ntuple_file_map() const
      { return ntuple_file_map_; }

  private:

    inline FilePropertiesManager() {
      this->load_file_properties();
    }

    inline void load_file_properties() {

      const char* path = std::getenv( "STV_ANALYSIS_DIR" );
      if ( path == nullptr ) throw std::runtime_error( "The environment"
        " variable STV_ANALYSIS_DIR is not set. Please set it and try again." );

      std::string in_file_name( path );
      in_file_name += "/file_properties.txt";
      std::ifstream in_file( in_file_name );

      std::string temp_line;
      while ( std::getline(in_file, temp_line) ) {
        // Ignore lines that begin with the '#' character (this allows for
        // comments in the normalization table file
        if ( temp_line.front() == '#' ) continue;

        // Read in the ntuple file name, the run number, and the file type from
        // the current line of the table file
        std::string file_name;
        int run;
        std::string type_str;
        std::istringstream temp_ss( temp_line );
        temp_ss >> file_name >> run >> type_str;

        // Convert the type string into an enum class value
        NtupleFileType type = string_to_file_type_map_.at( type_str );

        // If there is not an inner map for the current run number, then create
        // one
        if ( !ntuple_file_map_.count(run) ) {
          ntuple_file_map_[ run ] = std::map< NtupleFileType,
            std::set<std::string> >();
        }
        auto& run_map = ntuple_file_map_.at( run );

        // If there is not a set of files for the current ntuple file type,
        // then create one
        if ( !run_map.count(type) ) {
          run_map[ type ] = std::set< std::string >();
        }
        auto& file_set = run_map.at( type );

        // Insert the current file name into the appropriate place in the map
        file_set.insert( file_name );

        // For data files, also read in the trigger count and POT exposure
        // needed for normalization purposes
        if ( type == NtupleFileType::kOnBNB
          || type == NtupleFileType::kExtBNB )
        {
          int trigs;
          double pot;
          temp_ss >> trigs >> pot;

          // Store this information in the normalization map
          data_norm_map_[ file_name ] = TriggersAndPOT( trigs, pot );
        }
      }
    }

    // Outer keys are run numbers, inner keys are ntuple file types, values are
    // sets of file names
    std::map< int, std::map<NtupleFileType, std::set<std::string> > > ntuple_file_map_;

    // Keys are file names for processed data ntuples, values are objects
    // storing the corresponding number of triggers and POT exposure
    std::map< std::string, TriggersAndPOT > data_norm_map_;

    std::map< std::string, NtupleFileType > string_to_file_type_map_ = {
      { "onBNB", NtupleFileType::kOnBNB },
      { "extBNB", NtupleFileType::kExtBNB },
      { "numuMC", NtupleFileType::kNumuMC },
      { "nueMC", NtupleFileType::kIntrinsicNueMC },
      { "dirtMC", NtupleFileType::kDirtMC },
      { "detVarCV", NtupleFileType::kDetVarMCCV },
      { "detVarLYatten", NtupleFileType::kDetVarMCLYatten },
      { "detVarLYdown", NtupleFileType::kDetVarMCLYdown },
      { "detVarLYrayl", NtupleFileType::kDetVarMCLYrayl },
      { "detVarRecomb2", NtupleFileType::kDetVarMCRecomb2 },
      { "detVarSCE", NtupleFileType::kDetVarMCSCE },
      { "detVarWMAngleXZ", NtupleFileType::kDetVarMCWMAngleXZ },
      { "detVarWMAngleYZ", NtupleFileType::kDetVarMCWMAngleYZ },
      { "detVarWMdEdx", NtupleFileType::kDetVarMCWMdEdx },
      { "detVarWMX", NtupleFileType::kDetVarMCWMX },
      { "detVarWMYZ", NtupleFileType::kDetVarMCWMYZ },
    };

};