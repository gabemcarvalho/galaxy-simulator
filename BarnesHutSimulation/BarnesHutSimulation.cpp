#include "Globals.h"
#include "Config.h"
#include "SimulateBase.h"
#include "SimulateToyStar.h"
#include "SimulateSellwoodGalaxy.h"
#include "SimulateSellwoodSPHGalaxy.h"

int main()
{
    LoadFilenames(g_sFilenameFile);
    LoadConfig(g_sConfigFile.c_str());

    srand(g_iSeed);

    if (g_bUseToyStarModel)
    {
        RunToyStarSimulation();
    }
    else if (g_bUseSellwoodModel)
    {
        RunSellwoodGalaxySimulation();
    }
    else if (g_bUseSellwoodSPHModel)
    {
        RunSellwoodSPHGalaxySimulation();
    }
    else
    {
        RunBaseSimulation();
    }

    return 0;
}

