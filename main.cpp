#include <iostream>
#include <queue>
#include <string>

#include "Pendulum_Analysis.h"
#include "C:/Fitures/matrix/Matrix.h"
#include "msc.h"

const std::string file_extension = ".txt";

void smth(const std::string DataPath, std::string file_name, const size_t counter)
{
    Pendulum Pnd_arr[counter];
    for(size_t i = 1; i <= counter; ++i){
        Pnd_arr[i].Load_Pendulum_angle_history(DataPath + std::to_string(i) + "_r" + file_extension);
        Pnd_arr[i].Load_Pendulum_angle_history(DataPath + std::to_string(i) + file_extension);
            Pnd_arr[i].shift_angle_history_to_zero();
        Pnd_arr[i].find_freq_through_data(1024*2, 5);
        Pnd_arr[i].Save_freq(DataPath + std::to_string(i) + "_freq" + file_extension);
        Pnd_arr[i].Save_period(DataPath + std::to_string(i) + "_period" + file_extension);

        std::cout << i << "\n";
    }
    for(size_t i = 0; i < counter; ++i)
        Pnd_arr[i].info();
}

int main()
{
    std::string DataPath;
    std::string input_file_name;
    std::string output_file_name;

    DataPath = "";
    input_file_name = "01_r";

    smth("C:/Base/Inertia/RC_small/09.10.2020/Processing/", "", 29);

    /**Pendulum P1;
        P1.Load_Pendulum_angle_history(DataPath + input_file_name + file_extension);
            P1.shift_angle_history_to_zero();
    P1.find_freq_through_data(1024*2, 5);
    P1.Save_period(DataPath + input_file_name + "_period" + file_extension);*/

    //P1.info();
}
