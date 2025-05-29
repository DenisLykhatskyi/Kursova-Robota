#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define NUM_VARIANTS 3
#define INPUT_FILENAME "input.txt"
#define OUTPUT_FILENAME "output.txt"

const double PI = M_PI;

typedef struct {
    double T_total_time;
    double dt_time_step;
    double rho_density;
    double mu_viscosity;
    double A_const;
    double v0_initial_velocity;
    double d_diameter;
} VariantParameters;

double calculate_S_projected_area(double d_diameter) {
    return PI * d_diameter * d_diameter / 4.0;
}

double calculate_v_velocity(double t_current_time, double v0_initial_velocity, double A_const, double T_total_time) {
    double T_div_4 = T_total_time / 4.0;
    double Three_T_div_4 = 3.0 * T_total_time / 4.0;

    if (t_current_time >= 0 && t_current_time <= T_div_4) {
        return v0_initial_velocity + (4.0 * A_const / T_total_time) * t_current_time;
    } else if (t_current_time > T_div_4 && t_current_time <= Three_T_div_4) {
        return v0_initial_velocity + A_const;
    } else if (t_current_time > Three_T_div_4 && t_current_time <= T_total_time) {
        return v0_initial_velocity + A_const - (t_current_time - Three_T_div_4) * (4.0 * A_const / T_total_time);
    }
    return v0_initial_velocity;
}

double calculate_Re_reynolds_number(double v_velocity, double d_diameter, double rho_density, double mu_viscosity) {
    if (mu_viscosity == 0) return 0;
    return (v_velocity * d_diameter * rho_density) / mu_viscosity;
}

double calculate_phi_drag_coefficient(double Re_reynolds_number) {
    if (Re_reynolds_number <= 0) {
        return 0;
    }
    if (Re_reynolds_number <= 2.0) {
        return 24.0 / Re_reynolds_number;
    } else if (Re_reynolds_number <= 500.0) {
        return 18.5 / pow(Re_reynolds_number, 0.6);
    } else if (Re_reynolds_number <= 200000.0) { // 2*10^5
        return 0.44;
    } else {
        return 0.44; 
    }
}

double calculate_F_drag_force(double phi_drag_coefficient, double rho_density, double S_projected_area, double v_velocity) {
    return phi_drag_coefficient * rho_density * S_projected_area * v_velocity * v_velocity / 2.0;
}

int main() {
    VariantParameters variants[NUM_VARIANTS];
    FILE *inputFile, *outputFile;
    int i;
    double current_t;

    inputFile = fopen(INPUT_FILENAME, "r");
    if (inputFile == NULL) {
        perror("Помилка відкриття вхідного файлу");
        return 1;
    }

    for (i = 0; i < NUM_VARIANTS; ++i) {
        if (fscanf(inputFile, "%lf %lf %lf %lf %lf %lf %lf",
                   &variants[i].T_total_time, &variants[i].dt_time_step,
                   &variants[i].rho_density, &variants[i].mu_viscosity,
                   &variants[i].A_const, &variants[i].v0_initial_velocity,
                   &variants[i].d_diameter) != 7) {
            fprintf(stderr, "Помилка читання даних для варіанту %d з вхідного файлу.\n", i + 1);
            fclose(inputFile);
            return 1;
        }
    }
    fclose(inputFile);

    outputFile = fopen(OUTPUT_FILENAME, "w");
    if (outputFile == NULL) {
        perror("Помилка відкриття вихідного файлу");
        return 1;
    }

    for (i = 0; i < NUM_VARIANTS; ++i) {
        fprintf(outputFile, "Варіант %d:\n", i + 1);
        fprintf(outputFile, "Параметри: T=%.2f, dt=%.2f, rho=%.2f, mu=%.6e, A=%.2f, v0=%.2f, d=%.4f\n",
                variants[i].T_total_time, variants[i].dt_time_step, variants[i].rho_density,
                variants[i].mu_viscosity, variants[i].A_const, variants[i].v0_initial_velocity,
                variants[i].d_diameter);
        fprintf(outputFile, "Час (с)  Сила опору F (Н)\n");
        fprintf(outputFile, "--------------------------\n");

        VariantParameters current_params = variants[i];
        double S_area = calculate_S_projected_area(current_params.d_diameter);

        for (current_t = 0.0; current_t <= current_params.T_total_time + current_params.dt_time_step * 0.5; current_t += current_params.dt_time_step) {
            if (current_t > current_params.T_total_time && (current_t - current_params.T_total_time > current_params.dt_time_step * 0.1)) {
                break; 
            }
            if (current_t > current_params.T_total_time) {
                current_t = current_params.T_total_time;
            }

            double v = calculate_v_velocity(current_t, current_params.v0_initial_velocity, current_params.A_const, current_params.T_total_time);
            double Re = calculate_Re_reynolds_number(v, current_params.d_diameter, current_params.rho_density, current_params.mu_viscosity);
            double phi = calculate_phi_drag_coefficient(Re);
            double F = calculate_F_drag_force(phi, current_params.rho_density, S_area, v);

            fprintf(outputFile, "%-9.2f %.6f\n", current_t, F);
            

            if (fabs(current_t - current_params.T_total_time) < current_params.dt_time_step * 0.01) {
                break;
            }
        }
        fprintf(outputFile, "\n");
    }

    fclose(outputFile);
    printf("Обчислення завершено. Результати збережено у файл %s\n", OUTPUT_FILENAME);

    return 0;
}