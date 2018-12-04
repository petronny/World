//
// Created by muwang on 2017/11/24.
// extractsave input.wav f0_name.f0 sp_name.sp bap_name.bap frame_period fft_size f0_floor f0_ceil
// input.wav   : Input wave file
// frame_period: default: 5.0
// f0_floor    : default: 60.0
// f0_ceil     : default: 240.0
//

#include <stdio.h>
#include <stdlib.h>

// For .wav input/output functions.
#include "audioio.h"

// World core functions.
#include "world/d4c.h"
#include "world/dio.h"
#include "world/matlabfunctions.h"
#include "world/cheaptrick.h"
#include "world/stonemask.h"
#include "world/codec.h"

typedef struct {
    double frame_period;
    int fs;

    double *f0;
    double *time_axis;
    double f0_floor;
    double f0_ceil;
    int f0_length;

    double **spectrogram;
    double **bap;
    int fft_size;
} WorldParameters;

namespace {

    void F0EstimationDio(double *x, int x_length,
                         WorldParameters *world_parameters) {
        DioOption option = {0};
        InitializeDioOption(&option);

        // Modification of the option
        option.frame_period = world_parameters->frame_period;

        // Valuable option.speed represents the ratio for downsampling.
        // The signal is downsampled to fs / speed Hz.
        // If you want to obtain the accurate result, speed should be set to 1.
        option.speed = 1;

        // You can set the f0_floor below world::kFloorF0.
        option.f0_floor = world_parameters->f0_floor;
        option.f0_ceil  = world_parameters->f0_ceil;

        // You can give a positive real number as the threshold.
        // Most strict value is 0, but almost all results are counted as unvoiced.
        // The value from 0.02 to 0.2 would be reasonable.
        option.allowed_range = 0.1;

        // Parameters setting and memory allocation.
        world_parameters->f0_length = GetSamplesForDIO(world_parameters->fs,
                                                       x_length, world_parameters->frame_period);
        world_parameters->f0 = new double[world_parameters->f0_length];
        world_parameters->time_axis = new double[world_parameters->f0_length];
        double *refined_f0 = new double[world_parameters->f0_length];

        // printf("\nAnalysis\n");
        // DWORD elapsed_time = timeGetTime();
        Dio(x, x_length, world_parameters->fs, &option, world_parameters->time_axis,
            world_parameters->f0);
        // printf("DIO: %d [msec]\n", timeGetTime() - elapsed_time);

        // StoneMask is carried out to improve the estimation performance.
        // elapsed_time = timeGetTime();
        StoneMask(x, x_length, world_parameters->fs, world_parameters->time_axis,
                  world_parameters->f0, world_parameters->f0_length, refined_f0);
        // printf("StoneMask: %d [msec]\n", timeGetTime() - elapsed_time);

        for (int i = 0; i < world_parameters->f0_length; ++i)
            world_parameters->f0[i] = refined_f0[i];

        delete[] refined_f0;
    }

    void SpectralEnvelopeEstimation(double *x, int x_length,
                                    WorldParameters *world_parameters) {
        CheapTrickOption option = {0};
        // Note (2017/01/02): fs is added as an argument.
        InitializeCheapTrickOption(world_parameters->fs, &option);

        // Default value was modified to -0.15.
        option.q1 = -0.15;

        // Important notice (2017/01/02)
        // You can set the fft_size.
        // Default is GetFFTSizeForCheapTrick(world_parameters->fs, &option);
        // When fft_size changes from default value,
        // a replaced f0_floor will be used in CheapTrick().
        // The lowest F0 that WORLD can work as expected is determined
        // by the following : 3.0 * fs / fft_size.
        // option.f0_floor = 71.0;
        // option.fft_size = GetFFTSizeForCheapTrick(world_parameters->fs, &option);
        // We can directly set fft_size.
        option.fft_size = world_parameters->fft_size;

        // Parameters setting and memory allocation.
        world_parameters->spectrogram = new double *[world_parameters->f0_length];
        for (int i = 0; i < world_parameters->f0_length; ++i)
            world_parameters->spectrogram[i] =
                    new double[world_parameters->fft_size / 2 + 1];

        // DWORD elapsed_time = timeGetTime();
        CheapTrick(x, x_length, world_parameters->fs, world_parameters->time_axis,
                   world_parameters->f0, world_parameters->f0_length, &option,
                   world_parameters->spectrogram);
        // printf("CheapTrick: %d [msec]\n", timeGetTime() - elapsed_time);
    }

    void CodedAperiodicityEstimation(double *x, int x_length,
                                WorldParameters *world_parameters) {
        D4COption option = {0};
        InitializeD4COption(&option);

        // This parameter is used to determine the aperiodicity at 0 Hz.
        // If you want to use the conventional D4C, please set the threshold to 0.0.
        // Unvoiced section is counted by using this parameter.
        // Aperiodicity indicates high value when the frame is the unvoiced section.
        option.threshold = 0.85;

        // Parameters setting and memory allocation.
        double **aperiodicity = new double *[world_parameters->f0_length];
        for (int i = 0; i < world_parameters->f0_length; ++i)
            aperiodicity[i] = new double[world_parameters->fft_size / 2 + 1];

        // DWORD elapsed_time = timeGetTime();
        D4C(x, x_length, world_parameters->fs, world_parameters->time_axis,
            world_parameters->f0, world_parameters->f0_length,
            world_parameters->fft_size, &option, aperiodicity);
        // printf("D4C: %d [msec]\n", timeGetTime() - elapsed_time);

        // Code Aperiodicity.
        int number_of_aperiodicities = GetNumberOfAperiodicities(world_parameters->fs);
        world_parameters->bap = new double *[world_parameters->f0_length];
        for (int i = 0; i < world_parameters->f0_length; ++i)
            world_parameters->bap[i] = new double[number_of_aperiodicities];

	    CodeAperiodicity(aperiodicity, world_parameters->f0_length, world_parameters->fs,
                         world_parameters->fft_size, world_parameters->bap);

        for (int i = 0; i < world_parameters->f0_length; ++i)
            delete[] aperiodicity[i];
    }

    void FeatureWrite(const char *f0_path, const char *sp_path, const char *bap_path, WorldParameters *worldParameters){
        // Write f0
        FILE *f0_fp = fopen(f0_path, "wb");
        if (NULL == f0_fp) {
            printf("File cannot be opened: %s\n", f0_path);
            return;
        }
        fwrite(worldParameters->f0, 8, worldParameters->f0_length, f0_fp);
        fclose(f0_fp);

        int number_of_dimensions = worldParameters->fft_size / 2 + 1;
       // Write sp
        FILE *sp_fp = fopen(sp_path, "wb");
        if (NULL == sp_fp) {
            printf("File cannot be opened: %s\n", sp_path);
            return;
        }
        for (int i = 0; i < worldParameters->f0_length; i++){
            fwrite(worldParameters->spectrogram[i], 8, number_of_dimensions, sp_fp);
        }
        fclose(sp_fp);

        // Write bap
        FILE *bap_fp = fopen(bap_path, "wb");
        if (NULL == bap_fp) {
            printf("File cannot be opened: %s\n", bap_path);
            return;
        }
        int number_of_aperiodicities = GetNumberOfAperiodicities(worldParameters->fs);
        for (int i = 0; i < worldParameters->f0_length; i++){
            fwrite(worldParameters->bap[i], 8, number_of_aperiodicities, bap_fp);
        }
        fclose(bap_fp);
    }

}

int main(int argc, char *argv[]){
    if (argc != 9){
        printf("Error!\n");
        printf("Usage:\n");
        printf("\textractsave input.wav f0_name.f0 sp_name.sp bap_name.bap frame_period fft_size f0_floor f0_ceil\n");
        printf("Args:\n");
        printf("\t\tinput.wav:\tInput wave file path\n");
        printf("\t\tf0_name.f0:\tf0 save path\n");
        printf("\t\tsp_name.sp:\tsp save path\n");
        printf("\t\tbap_name.ap:\tbap save path\n");
        printf("\t\tframe_period:\tframe shift [ms]\n");
        printf("\t\tfft_size:\tfft size [2^N]\n");
        printf("\t\tf0_floor:\tf0 floor [double]\n");
        printf("\t\tf0_ceil:\tf0 ceil [double]\n");
        return -2;
    }

    // Memory allocation is carried out in advanse.
    // This is for compatibility with C language.
    int x_length = GetAudioLength(argv[1]);
    if (x_length <= 0) {
        if (x_length == 0) printf("error: File not found.\n");
        else printf("error: The file is not .wav format.\n");
        return -1;
    }
    double *x = new double[x_length];
    // wavread() must be called after GetAudioLength().
    int fs, nbit;
    wavread(argv[1], &fs, &nbit, x);
    // DisplayInformation(fs, nbit, x_length);

    //---------------------------------------------------------------------------
    // Analysis part
    //---------------------------------------------------------------------------
    WorldParameters world_parameters = { 0 };
    // You must set fs and frame_period before analysis/synthesis.
    world_parameters.fs = fs;
    world_parameters.frame_period = atof(argv[5]);
    world_parameters.fft_size = atoi(argv[6]);
    world_parameters.f0_floor = atof(argv[7]);
    world_parameters.f0_ceil = atof(argv[8]);

    // F0 estimation
    // DIO
    F0EstimationDio(x, x_length, &world_parameters);

    // Spectral envelope estimation
    SpectralEnvelopeEstimation(x, x_length, &world_parameters);

    // Coded aperiodicity estimation by D4C
    CodedAperiodicityEstimation(x, x_length, &world_parameters);

    // Write feature
    FeatureWrite(argv[2], argv[3], argv[4], &world_parameters);

    return 0;
}