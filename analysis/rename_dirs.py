''' Convert old directory organization to new (Dec2014) format '''

import argparse
import os
import shutil
import glob



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='.')
    parser.add_argument('--name', type=str, required=True, help='Name of subdirectory.')
    args = parser.parse_args()

    name = args.name

    ## Tf_#/Ti_Tf_dT.txt       ->  iteration_#/short_Ti_Tf_dT
    ## Tf_#/Tf.txt             ->  iteration_#/short_Tf
    ## Tf_#/T_array.txt        ->  iteration_#/short_temps
    ## Tf_#/T_array_last.txt   ->  iteration_#/short_temps_last
    ## Tf_#/whamQ              ->  iteration_#/short_wham

    ## Mut_#/Ti_Tf_dT.txt      ->  iteration_#/long_Ti_Tf_dT
    ## Mut_#/Tf.txt            ->  iteration_#/long_Tf
    ## Mut_#/T_array.txt       ->  iteration_#/long_temps
    ## Mut_#/T_array_last.txt  ->  iteration_#/long_temps_last
    ## Mut_#/whamQ             ->  iteration_#/long_wham

    cwd = os.getcwd()
    os.chdir(name)
    
    print "Moving directory:"
    for i in range(10):
        y = "iteration_%d" % i
        x1 = "Tf_%d" % i
        if os.path.exists(x1):
            print "  %s" % x1
            if not os.path.exists(y):
                os.mkdir(y)

            if os.path.exists("%s/Ti_Tf_dT.txt" % x1):
                shutil.move("%s/Ti_Tf_dT.txt" % x1,"%s/short_Ti_Tf_dT" % y)
            if os.path.exists("%s/Tf.txt" % x1):
                shutil.move("%s/Tf.txt" % x1,"%s/short_Tf" % y)
            if os.path.exists("%s/T_array.txt" % x1):
                shutil.move("%s/T_array.txt" % x1,"%s/short_temps" % y)
            if os.path.exists("%s/T_array_last.txt" % x1):
                shutil.move("%s/T_array_last.txt" % x1,"%s/short_temps_last" % y)
            if os.path.exists("%s/whamQ" % x1):
                shutil.move("%s/whamQ" % x1,"%s/short_wham" % y)

            other_data =  glob.glob("%s/*" % x1)
            for data in other_data:
                shutil.move(data,"%s/" % y)

        x2 = "Mut_%d" % i
        if os.path.exists(x2):
            print "  %s" % x2
            if not os.path.exists(y):
                os.mkdir(y)

            if os.path.exists("%s/Ti_Tf_dT.txt" % x2):
                shutil.move("%s/Ti_Tf_dT.txt" % x2,"%s/long_Ti_Tf_dT" % y)
            if os.path.exists("%s/Tf.txt" % x2):
                shutil.move("%s/Tf.txt" % x2,"%s/long_Tf" % y)
            if os.path.exists("%s/T_array.txt" % x2):
                shutil.move("%s/T_array.txt" % x2,"%s/long_temps" % y)
            if os.path.exists("%s/T_array_last.txt" % x2):
                shutil.move("%s/T_array_last.txt" % x2,"%s/long_temps_last" % y)
            if os.path.exists("%s/whamQ" % x2):
                shutil.move("%s/whamQ" % x2,"%s/long_wham" % y)

            other_data =  glob.glob("%s/*" % x2)
            for data in other_data:
                shutil.move(data,"%s/" % y)

    os.chdir(cwd)
