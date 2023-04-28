#ifndef INCLUDE_OSCAROUTPUT_H_
#define INCLUDE_OSCAROUTPUT_H_

#include "smash/file.h"
#include "smash/oscaroutput.h"
#include "smash/outputparameters.h"
#include "smash/particles.h"

#include "gen.h"
#include "params.h"

void write_oscar_output() {
  // initialize Oscar output
  const std::filesystem::path OutputPath = params::sSpectraDir;
  std::unique_ptr<smash::OutputInterface> OscarOutput =
       create_oscar_output("Oscar2013", "Particles", OutputPath,
                           smash::OutputParameters());

  for(int iev=0; iev<params::NEVENTS; iev++){
    // Create Particles oject which contains the full particle list of each event
    std::unique_ptr<smash::Particles> particles = std::make_unique<smash::Particles>();

    for(int inpart=0; inpart<gen::npart[iev]; inpart++) {
      smash::ParticleData* data = gen::pList[iev][inpart];
      particles->insert(*data);
    }
    // Set dummy values for event info:
    // {impact parameter, box_length, current time, total_kinetic_energy, total mean field energy, number of test particles, empty event (no collisions?)}
    smash::EventInfo event_info{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1, false};
    // Write output
    OscarOutput->at_eventend(*particles, iev, event_info);
  }
}

#endif // INCLUDE_OSCAROUTPUT_H_
