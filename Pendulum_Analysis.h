#ifndef PENDULUM_ANALYSIS_H_INCLUDED
#define PENDULUM_ANALYSIS_H_INCLUDED

#include "C:/Fitures/matrix/Matrix.h"
#include <string>

class Pendulum
{
private:
protected:
public:

    enum _MODE_PEAK_FREQ
    {
        SIMPLE,
        POLYNOM_APPROX2
    };

    enum _MODE_LOAD_FILE
    {
        PROP_FILE,
        LIR_PROP_FILE,
        LIR_ROW_FILE
    };

    ///input
    double* angle_history;
    size_t data_length;
    double discr_t;
    ///output
    double discr_freq;
    uint16_t zeroes_fitting_factor;
    double* freq;
    size_t freq_length;
    ///information
    std::string id;///assosiation with inertia Oscillation class etc
    const std::string file_extension;
    ///MODES
    _MODE_PEAK_FREQ MODE_PEAK_FREQ = _MODE_PEAK_FREQ::SIMPLE;
    _MODE_LOAD_FILE MODE_LOAD_FILE = _MODE_LOAD_FILE::LIR_PROP_FILE;


    Pendulum();
    Pendulum(double* n_angle_history, const size_t n_angle_history_length, const double &n_discr_t);

    void set_angle_history(double* n_angle_history, const size_t n_angle_history_length, const double &n_discr_t);
    void set_freq(double* n_freq, const size_t n_freq_length, const double &n_discr_freq);
    void set_id(const std::string &n_id);

    double count_window_freq(size_t index_from, size_t window_length);
    void find_freq_through_data(size_t window_size, size_t window_step);

    size_t length_of_angle_history_file(std::string);
    void fix_LIR_ofscale(double *angle_history, size_t data_length);
    void Load_Pendulum_angle_history(const std::string file_name);
    void shift_angle_history_to_zero();

    void find_pendulum_frequency();
    void Correct_Pendulum_Frequency();/// using polynomial correction depends on amplitude

    void info();
};


void Load_Data(const double* data, const size_t data_length);

void Correct_Pendulum_Frequency(const double* data, const size_t data_length);

void find_max_spectrum_element(double* spectrum, const size_t data_length, double &max_value, size_t &max_value_index);


//void manage

#endif // PENDULUM_ANALYSIS_H_INCLUDED
