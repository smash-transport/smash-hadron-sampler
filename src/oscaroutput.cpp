#include "smash/oscaroutput.h"
#include "smash/file.h"
#include "smash/outputparameters.h"
#include "smash/particles.h"

#include "gen.h"
#include "oscaroutput.h"
#include "params.h"

void write_oscar_output() {
  // initialize Oscar output
  const std::filesystem::path OutputPath = params::output_directory;
  std::unique_ptr<smash::OutputInterface> OscarOutput = create_oscar_output(
      "Oscar2013", "Particles", OutputPath, smash::OutputParameters());

  for (int iev = 0; iev < params::number_of_events; iev++) {
    // Create Particles oject which contains the full particle list of each
    // event
    std::unique_ptr<smash::Particles> particles =
        std::make_unique<smash::Particles>();

    for (int inpart = 0; inpart < gen::npart[iev]; inpart++) {
      smash::ParticleData *data = gen::pList[iev][inpart];
      particles->insert(*data);
    }
    // Set dummy values for event info:
    // {impact parameter, box_length, current_time, E_kinetic_total, total
    // mean field energy, E_total, number of test particles, n_ensembles,
    // !projectile_target_interact, kinematic_cut_for_SMASH_IC}
    // See src/experiment.cc in the smash repository for more detailed
    // information on event_info
    smash::EventInfo event_info{0.0, 0.0, 0.0, 0.0,   0.0,
                                0.0, 1,   0,   false, false};
    smash::EventLabel event_id{iev, 0};
    // Write output
    OscarOutput->at_eventend(*particles, event_id, event_info);
  }
}
