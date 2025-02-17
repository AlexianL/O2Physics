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
#include <vector>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"

#include "MathUtils/detail/TypeTruncation.h"

#include "PWGCF/DataModel/CorrelationsDerived.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::math_utils::detail;

// #define FLOAT_PRECISION 0xFFFFFFF0
#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct FilterMft {
  O2_DEFINE_CONFIGURABLE(cfgVerbosity, int, 0, "Verbosity level (0 = major, 1 = per collision)")

  //  configurables for MFT tracks
  Configurable<double> etaMftTrackMax{"cfgEtaMftTrackMax", -2.4, "Maximum value for the eta of MFT tracks"};
  Configurable<double> etaMftTrackMin{"cfgEtaMftTrackMin", -5, "Minimum value for the eta of MFT tracks"};
  Configurable<int> nClustersMftTrack{"cfgNClustersMftTrack", 5, "Minimum number of clusters for the reconstruction of MFT tracks"};

  Produces<aod::CFMftTracks> outputMftTracks;

  template <typename TTrack>
  bool isAcceptedMftTrack(TTrack const& mftTrack)
  {
    // cut on the eta of MFT tracks
    if (mftTrack.eta() > etaMftTrackMax || mftTrack.eta() < etaMftTrackMin) {
      return false;
    }

    // cut on the number of clusters of the reconstructed MFT track
    if (mftTrack.nClusters() < nClustersMftTrack) {
      return false;
    }

    return true;
  }

  void processData(aod::Collisions::iterator const&, aod::BCsWithTimestamps const&, aod::CFCollRefs const& cfcollisions, aod::CFTrackRefs const& cftracks, aod::MFTTracks const& mfttracks)
  {
    if (cfcollisions.size() <= 0 || cftracks.size() <= 0) {
      return; // rejected collision
    }

    if (cfgVerbosity > 0 && mfttracks.size() > 0) {
      LOGF(info, "MFT Tracks for collision: %lu, cfcollisions: %lu, CFTracks: %lu", candidates.size(), cfcollisions.size(), cftracks.size());
    }

    for (auto& mfttrack : mfttracks) {

      if (isAcceptedMftTrack(mfttrack)) {
        outputMftTracks(cfcollisions.begin().globalIndex(),
                        mfttrack.pt(), mfttrack.eta(), mfttrack.phi(), mfttrack.nClusters());
      }
    }
  }
  PROCESS_SWITCH(FilterMft, processData, "Process data MFT tracks", true);
}; // struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FilterMft>(cfgc, TaskName{"filter-mft"})};
}
