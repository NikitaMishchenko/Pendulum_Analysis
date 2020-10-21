#include <iostream>
#include <fstream>
#include <algorithm>
#include <math.h>

#include "Pendulum_Analysis.h"
#include "C:/Fitures/FFT/FFT.h"
#include "C:/Fitures/matrix/Matrix.h"

Pendulum::Pendulum()
{
    angle_history = nullptr;
    data_length = 0;
    discr_t = 0.0;

    discr_freq = 0.0;
    zeroes_fitting_factor = 0;
    freq = nullptr;
    freq_length = 0;

    id = "0000";
}

Pendulum::Pendulum(double* n_angle_history, const size_t n_angle_history_length, const double &n_discr_t)
{
    angle_history = new double [n_angle_history_length];
    for(size_t i = 0; i < n_angle_history_length; ++i){
        angle_history[i] = n_angle_history[i];
    }
    data_length = n_angle_history_length;
    discr_t = n_discr_t;

    discr_freq = 0.0;
    zeroes_fitting_factor = 0;
    freq = nullptr;
    freq_length = 0;

    id = "0000";
}

void Pendulum::set_angle_history(double* n_angle_history, const size_t n_angle_history_length, const double &n_discr_t)
{
    delete [] angle_history;

    angle_history = new double [n_angle_history_length];
    for(size_t i = 0; i < n_angle_history_length; ++i){
        angle_history[i] = n_angle_history[i];
    }
    data_length = n_angle_history_length;
    discr_t = n_discr_t;
}

void Pendulum::set_freq(double* n_freq, const size_t n_freq_length, const double &n_discr_freq)
{
    delete [] freq;

    freq = new double [n_freq_length];
    for(size_t i = 0; i < n_freq_length; ++i){
        freq[i] = n_freq[i];
    }
    freq_length = n_freq_length;
    discr_freq = n_discr_freq;
}

void Pendulum::set_id(const std::string &n_id)
{
    id = n_id;
}

double Pendulum::count_window_freq(size_t index_from, size_t window_length){
    double* window = new double [window_length];
    double* w_spectrum = new double [window_length];

    for(size_t i = index_from; i < index_from+window_length; ++i){
        window[i-index_from] = angle_history[i];
        w_spectrum[i-index_from] = 0.0;
    }

    int zeroes_fitting_factor = 2;
    ///depends on enum mode

    static uint16_t number_of_w = 0;
    ///static?///
    simple_FFT sf(window, window_length, discr_t, zeroes_fitting_factor);
        sf.general_FFT();

    std::ofstream fout_s("window_" + std::to_string(number_of_w) + ".txt");
        for(size_t i = 0; i < window_length; ++i){
            fout_s << i*sf.discr_t << "\t" << sf.signal[i] << std::endl;
        }
    fout_s.close();

    std::ofstream fout_w("window_fft_" + std::to_string(number_of_w) + ".txt");
        for(size_t i = 0; i < window_length; ++i){
            fout_w << i/sf.discr_t/sf.Nft << "\t" << log(sf.power[i]) << std::endl;
        }
    fout_w.close();
    number_of_w++;

    double peak_freq = 0.0;
        if(MODE_PEAK_FREQ == _MODE_PEAK_FREQ::SIMPLE){
            size_t max_element_number = std::distance (sf.power, std::max_element(sf.power, sf.power+window_length));
            //std::cout << "max_spectrum_element = " << max_element_number << "\tpeak freq = " << 1.0/sf.discr_t/sf.Nft*max_element_number << std::endl;
            peak_freq = 1.0/sf.discr_t/sf.Nft*max_element_number;
            return peak_freq;
        }
        if(MODE_PEAK_FREQ == _MODE_PEAK_FREQ::SIMPLE){
            ///approximate and get value
            peak_freq = 0.0;
            return peak_freq;
        }

    delete [] window;
    delete [] w_spectrum;
    return peak_freq;
}

void Pendulum::find_freq_through_data(size_t window_size, size_t window_step){
    size_t n_freq_length = data_length/window_size + 1;
    delete [] freq;
    freq = new double [n_freq_length];
    std::cerr << n_freq_length << std::endl;
    size_t k = 0;
    for(size_t i = 0; i < data_length - window_size; )
    {
        freq[k] = this->count_window_freq(i, window_size);
        i += window_step;
        k++;
    }
}

size_t Pendulum::length_of_angle_history_file(std::string file_name){
    std::ifstream fin(file_name);
    std::string buff = "";
    int height = 0;
    while(!fin.eof())
    {
        getline(fin, buff);
        height++;
    }
    fin.close();
    return height;
}

void Pendulum::fix_LIR_ofscale(double *angle_history, size_t data_length){
    double gap = 0.0, border = 350.0;
    bool flag = false; ///нужно ли выполнять смещение участка?
    for(size_t i = data_length-2; i > 0; --i)
    {
        gap = fabs(angle_history[i+1] - angle_history[i]);

        if(gap >= border)
            flag = true;

        if(gap < border)
            flag = false;

        if(flag==true)
            angle_history[i] += 360.0;
    }
}

void Pendulum::Load_Pendulum_angle_history(std::string file_name){
    size_t n_data_length = length_of_angle_history_file(file_name);
    if(data_length != n_data_length){

        data_length = n_data_length;
        delete angle_history;
        angle_history = new double [data_length];
    }

    std::ifstream fin(file_name);
    std::cerr << file_name << std::endl;
    double buff = 0.0;
    size_t counter = 0;
    fin >> buff;
        discr_t = buff;
    fin >> buff;
        angle_history[counter] = buff;
    counter++;
    fin >> buff;
        discr_t = buff - discr_t;
    fin >> buff;
        angle_history[counter] = buff;

    while(counter <= data_length){
        fin >> buff;
        fin >> buff;
        angle_history[counter] = buff;
        counter++;
    }

    if(MODE_LOAD_FILE == _MODE_LOAD_FILE::LIR_PROP_FILE)
    {
        fix_LIR_ofscale(angle_history, data_length);
        std::cerr << "DATA FIXED\n";
    }
}

void Pendulum::shift_angle_history_to_zero()///FIXME
{
    for(size_t i = 0; i < data_length; i++)
        angle_history[i] -= angle_history[data_length-1];
    //double b_a = angle_history[data_length-1];
    //std::transform(angle_history, angle_history+data_length, angle_history, [&b_a](double current){return (current - b_a);});
};

void Pendulum::info()
{
    std::cerr << "Pendulum object:\n";
    std::cerr << "\tid = " << id << std::endl;
    std::cerr << "\tdata_length = " << data_length << std::endl;
    std::cerr << "\tfreq_length = " << freq_length << std::endl;
}

void foo(double value)
{
    std::cout << value << "\t";
}




void find_max_spectrum_element(double* spectrum, const size_t data_length, double &max_value, size_t &max_value_index)
{
    std::for_each(&spectrum[0], &spectrum[0+data_length], &foo);
    //std::cout << std::endl;
    std::cout << *std::max_element(&spectrum[0],&spectrum[0+data_length]);

};

