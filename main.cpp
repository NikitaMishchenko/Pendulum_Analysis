#include <iostream>
#include <queue>
#include <string>

#include "Pendulum_Analysis.h"
#include "C:/Fitures/matrix/Matrix.h"
#include "msc.h"

const std::string file_extension = ".txt";

void proceed_list_pendulum_procedure(const std::string DataPath, std::string file_name, const size_t from, const size_t to)
{

    size_t counter = to - from;
    Pendulum Pnd_arr[counter];
    for(size_t i = 0; i <= counter; ++i){
        Pnd_arr[i].Load_Pendulum_angle_history(DataPath + std::to_string(i+from) + "_r" + file_extension);

        Pnd_arr->MODE_REPORT = Pnd_arr->_MODE_REPORT::EVERYTHING;
        //Pnd_arr[i].Load_Pendulum_angle_history(DataPath + std::to_string(i) + file_extension);
            Pnd_arr[i].shift_angle_history_to_zero();
        Pnd_arr[i].find_freq_through_data(1024*2, 5);

        ///SAVE
            Pnd_arr[i].Save_all_angle_history(DataPath + std::to_string(i+from) + "_all_angle_history" + file_extension);
            Pnd_arr[i].Save_freq(DataPath + std::to_string(i+from) + "_freq" + file_extension);
            Pnd_arr[i].Save_period(DataPath + std::to_string(i+from) + "_period" + file_extension);
            Pnd_arr[i].Save_envelop(DataPath + std::to_string(i+from) + "_envelop" + file_extension);
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

    //proceed_list_pendulum_procedure("C:/Base/Inertia/RC_small/09.10.2020/Processing/", "", 11, 28);

    /**Pendulum P1;
        P1.Load_Pendulum_angle_history(DataPath + input_file_name + file_extension);
            P1.shift_angle_history_to_zero();
    P1.find_freq_through_data(1024*2, 5);
    P1.Save_period(DataPath + input_file_name + "_period" + file_extension);*/

    //P1.info();
}
