#include <iostream>
#include <queue>
#include <string>

#include "Pendulum_Analysis.h"
#include "C:/Fitures/matrix/Matrix.h"

const std::string file_extension = ".txt";

void Load_Oscillation(const std::string DataPath, const std::string input_file_name, const std::string output_file_name)
{

}



int main()
{
    //size_t data_length = 11, max_value_index = 1;
    //double max_value = 0.0;
    //double spectrum[] = {1,2,3,4,5,7,8,9,7,5,0};
    //find_max_spectrum_element(spectrum, data_length, max_value, max_value_index);

    matrix Data;
    std::string DataPath;
    std::string input_file_name;
    std::string output_file_name;

    DataPath = "";
    input_file_name = "01";

    Pendulum P1;

    //Data.Load_matrix_emptstringend(DataPath + input_file_name + file_extension, 1);
    P1.Load_Pendulum_angle_history(DataPath + input_file_name + file_extension);

    P1.shift_angle_history_to_zero();
    ///P1.MODE_MAX_FREQ = Pendulum::_MODE_MAX_FREQ::POLYNOM_APPROX2;
    //std::cout << P1.MODE_PEAK_FREQ << std::endl;
    P1.find_freq_through_data(1024*8, 1000);
    //P1.count_window_freq(5000, 24000);

    //P1.info();
}
