#include "model_param_choice.h"

std::vector<bool> to_binary(int num_to_convert_to_binary, int num_bits_in_out_vec)
{
    std::vector<bool> r;

    // make binary vec of minimum size backwards (LSB at .end() and MSB at .begin())
    while (num_to_convert_to_binary > 0)
    {
        //cout << " top of loop" << endl;
        if (num_to_convert_to_binary % 2 == 0)
            r.push_back(0);
        else
            r.push_back(1);
        num_to_convert_to_binary = num_to_convert_to_binary / 2;
    }

    while(r.size() < num_bits_in_out_vec)
        r.push_back(0);

    std::reverse(r.begin(),r.end());

    return r;
}


std::vector<bool> current_model_choice(std::string  model_choice) {
    std::vector<std::vector<bool>> models(4);
    int cases = log2(models.size());

    for (int i = 0; i < models.size(); i++)
    {
        models[i] = to_binary(i, cases);
    }
    std::map<std::string ,std::vector<bool>> mymap;
    
    mymap["eq"] = models[0];
    mymap["r_dif"] = models[1];
    mymap["a_dif"] = models[2];
    mymap["dif"] = models[3];

    return mymap[model_choice];
}