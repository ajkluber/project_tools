import newsystem

def get_system(opts):
    System = newsystem.NewSystem(opts["Subdir"],Tf_it=opts["Tf_Iteration"],
                                Tf_act_dir=opts["Tf_Active_Directory"],Tf_refine=opts["Tf_Refinements"],
                                Mut_it=opts["Mut_Iteration"],Mut_act_dir=opts["Mut_Active_Directory"])
    return System

def new_system(subdir):
    System = newsystem.NewSystem(subdir)
    return System

def load_systems(subdirs):
    ''' Create systems from saved options in system.info'''
    Systems = []
    for subdir in subdirs:
        print "Loading model from subdirectory: ", subdir
        System = load_system(subdir)
        Systems.append(System)
    return Systems

def new_systems(subdirs):
    ''' Create systems from saved options in system.info'''
    Systems = []
    for subdir in subdirs:
        print "Loading model from subdirectory: ", subdir
        System = new_system(subdir)
        Systems.append(System)
    return Systems

def load_system(subdir):
    ''' Given subdir that contains model.info options file. Read in options and
        create corresponding model.'''
    info_file = open(subdir+'/system.info','r')
    line = info_file.readline()
    options = {}
    while line != '':
        print line
        field = line.split()[1]
        value = info_file.readline()
        if field == "Tf_refinements":
            temp = [int(value[:-1].split()[0])]
            line = info_file.readline()
            print line
            while (line[:1] != '[') and (line != ''):
               temp.append(int(line[:-1].split()[0]))
            options[field] = temp
        else:
            options[field] = value[:-1]
            line = info_file.readline()
    options = check_fields(options)
    System = get_system(options)
    return System

def check_fields(fields):
    ''' Check the integrity of the system information'''

    print "Checking system information..."
    print fields

    ## For backwards compatibility. Convert old names to new names.
    old_fields = { "Main_Path":"Path","subdir":"Subdir",
                   "Tf_iteration":"Tf_Iteration",
                   "Tf_active_directory":"Tf_Active_Directory",
                   "Tf_refinements":"Tf_Refinements",
                   "mutation_iteration":"Mut_Iteration",
                   "mutation_active_directory":"Mut_Active_Directory"}
    
    ## Replace old style keys with new style keys. 
    for key in fields.keys():
        if key in old_fields.keys():
            fields[old_fields[key]] = fields[key]

#    fields["Path"] = fields["Main_Path"]
#    fields["Subdir"] = fields["subdir"]
#    fields["Tf_Iteration"] = fields["Tf_iteration"]
#    fields["Tf_Active_Directory"] = fields["Tf_active_directory"]
#    fields["Tf_Refinements"] = fields["Tf_refinements"]
#    fields["Mut_Iteration"] = fields["mutation_iteration"]
#    fields["Mut_Active_Directory"] = fields["mutation_active_directory"]

    print "Using system information.."
    print fields

    return fields
