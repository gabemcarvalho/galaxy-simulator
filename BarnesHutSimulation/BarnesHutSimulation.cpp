#include "Globals.h"
#include "Config.h"
#include "SimulateBase.h"
#include "SimulateToyStar.h"

int main()
{
    LoadFilenames(g_sFilenameFile);
    LoadConfig(g_sConfigFile.c_str());

    srand(g_iSeed);

    if (g_bUseToyStarModel)
    {
        RunToyStarSimulation();
    }
    else
    {
        RunBaseSimulation();
    }

    return 0;
}

