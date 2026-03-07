#include "QGSP_BERT.hh"
#include "G4Gamma.hh"
#include "G4ProcessManager.hh"

class CaloXPhysicsList : public QGSP_BERT
{
public:
    CaloXPhysicsList() : QGSP_BERT() {}
    ~CaloXPhysicsList() {}

    void ConstructProcess() override
    {
        QGSP_BERT::ConstructProcess(); // Call base class to set up standard processes

        // Now, remove the photonNuclear process
        G4ParticleDefinition *gamma = G4Gamma::GammaDefinition();
        G4ProcessManager *pmanager = gamma->GetProcessManager();
        if (pmanager)
        {
            G4int nProcesses = pmanager->GetProcessListLength();
            for (G4int i = 0; i < nProcesses; i++)
            {
                G4VProcess *process = (*pmanager->GetProcessList())[i];
                if (process->GetProcessName() == "photonNuclear")
                {
                    pmanager->RemoveProcess(process);
                    G4cout << "Removed process: " << process->GetProcessName() << G4endl;
                    break;
                }
            }
        }
    }
};