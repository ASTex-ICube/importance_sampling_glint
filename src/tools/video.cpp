#include <fstream>
#include <string>

#include "geometry.h"
#include "pbrt.h"

using namespace pbrt;

// t : time variable
// T : total time
// x : acceleration time
Float acc(Float t, Float T, Float x) {
    if (t >= 0.f && t < x) return t * t * 0.5f / (T * x - x * x);
    if (t >= x && t < T - x) return (t * x - x * x * 0.5f) / (T * x - x * x);
    if (t >= T - x && t <= T + ShadowEpsilon)
        return (-t * t * 0.5f + T * t - T * T * 0.5f + T * x - x * x) /
               (T * x - x * x);
    return 0.f;
}

static void usage() {
    fprintf(stderr,
            "usage: video <input scene filename> <output filename>\n<input "
            "scene filename>\nnumber of variables varying over time\ninitial "
            "value for variable 0\nfinal value for variable 0\ninitial value "
            "for variable 1\nfinal value for variable 1\n...\n...\ntotal "
            "number of frames\nnumber of frames for "
            "acceleration/slowdown\nnumber of threads\npython3\n"
            "scene pbrt\n");
    exit(1);
}

int main(int argc, char *argv[]) {
    google::InitGoogleLogging(argv[0]);
    FLAGS_stderrthreshold = 1;  // Warning and above.

    std::string paramsFileName, outputFileName;
    for (int i = 1; i < argc; ++i) {
        if (!strcmp(argv[i], "--help") || !strcmp(argv[i], "-h"))
            usage();
        else if (paramsFileName.empty())
            paramsFileName = argv[i];
        else if (outputFileName.empty())
            outputFileName = argv[i];
        else
            usage();
    }

    std::ifstream paramsStream(paramsFileName);

    if (!paramsStream) {
        fprintf(stderr, "can't open %s", paramsFileName.c_str());
        return 1;
    }

    paramsStream.exceptions(paramsStream.exceptions() | std::ios_base::badbit);

    unsigned int nVariables;

    paramsStream >> nVariables;

    std::vector<Float> initialValues(nVariables);
    std::vector<Float> finalValues(nVariables);

    for (unsigned int i = 0; i < nVariables; ++i) {
        paramsStream >> initialValues[i];
        paramsStream >> finalValues[i];
    }

    unsigned int totalNFrames;
    unsigned int nFramesAcc;
    int nThreads;  // number of threads
    bool python3;
    std::string scene;

    paramsStream >> totalNFrames;
    paramsStream >> nFramesAcc;
    CHECK_GT(totalNFrames, nFramesAcc * 2);
    paramsStream >> nThreads;
    paramsStream >> python3;

    Float T = totalNFrames / 24.f;
    Float tAcc = nFramesAcc / 24.f;
    std::cout << "T: " << T << std::endl;
    std::cout << "tAcc: " << tAcc << std::endl;

    char cChar;
    while (paramsStream.get(cChar)) {
        scene += cChar;
    }

    std::ofstream launcher("launcher.py");
    launcher << "import subprocess" << std::endl;

    for (int frame = 0; frame < totalNFrames; ++frame) {
        Float t = frame / 24.f;

        std::ostringstream index;
        if (frame < 10)
            index << "_000" << frame;
        else if (frame < 100)
            index << "_00" << frame;
        else if (frame < 1000)
            index << "_0" << frame;
        else
            index << "_" << frame;

        Float accValue = acc(t, T, tAcc);

        std::string::size_type n;
        std::string::size_type i = 0;
        unsigned int j = 0u;

        std::string sceneFrameI(scene);

        n = sceneFrameI.find("%imageName", i);
        if (n != std::string::npos) {
            sceneFrameI.replace(n, 10, outputFileName + index.str());
        }

        do {
            n = sceneFrameI.find("%f", i);
            if (n != std::string::npos) {
                CHECK_LT(j, nVariables);

                Float initialValueJ = initialValues[j];
                Float finalValueJ = finalValues[j];
                Float intervalJ = finalValueJ - initialValueJ;
                std::ostringstream valueJString;
                valueJString << initialValueJ + accValue * intervalJ;
                sceneFrameI.replace(n, 2, valueJString.str());
            } else
                break;
            i = n + 1;
            ++j;
        } while (true);
        CHECK_EQ(j, nVariables);

        std::string sceneFileNameI(outputFileName + index.str() + ".pbrt");
        std::string lofFileNameI(outputFileName + index.str() + ".log");
        std::string scriptFileNameI(outputFileName + index.str() + ".py");

        std::ofstream oF(sceneFileNameI);
        if (!oF) {
            fprintf(stderr, "can't open output file %s",
                    sceneFileNameI.c_str());
            return 1;
        }

        oF << sceneFrameI;

        std::ofstream oFSh(scriptFileNameI);
        if (!oFSh) {
            fprintf(stderr, "can't open output file %s",
                    scriptFileNameI.c_str());
            return 1;
        }

        oFSh << "import subprocess" << std::endl;
        oFSh << "process = subprocess.run(['pbrt','--nthreads',str(" << nThreads
             << "),'" << sceneFileNameI
             << "'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, "
                "universal_newlines=True)"
             << std::endl;

        oFSh << "output_file = open(\"" << lofFileNameI << "\", \"w\")"
             << std::endl;
        oFSh << "output_file.write(process.stdout)" << std::endl;
        oFSh << "output_file.close()" << std::endl;

        launcher << "process = subprocess.run(['python";
        if(python3)
            launcher << "3";
        launcher << "','" << scriptFileNameI
                 << "'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, "
                    "universal_newlines=True)"
                 << std::endl;
        launcher << "print(process.stdout)" << std::endl;
    }

    return 0;
}
