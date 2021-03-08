import os

location = "src/SusyLesHouches.cc"
temp_location = location + "_temp"

lines = open(location, 'r').readlines()

# Find where the GAMBIT patch ends.
linenum = 0
with open(location) as f:
    for num, line in enumerate(f, 1):
        if "(blockName == \"nmssmrun\") ifail=nmssmrun.set(linestream)" in line:
            linenum = num+1
            break

with open(temp_location, 'w') as f:
    # Write the stuff at the beginning...
    for i in range(linenum):
        f.write(lines[i])
    # Write the source specific to the model...
    f.write((
        "      // LH blocks added by GUM\n"
        "      if (blockName == \"dmint\") ifail=dmint.set(linestream);\n"
        "\n"
    ))
    # Then write the rest. Voila: the cheap man's patch.
    for i in range(len(lines)-linenum):
        f.write(lines[i+linenum])

os.remove(location)
os.rename(temp_location, location)

location = "include/Pythia8/SusyLesHouches.h"
temp_location = location + "_temp"

lines = open(location, 'r').readlines()

# Find where the GAMBIT patch ends.
linenum = 0
with open(location) as f:
    for num, line in enumerate(f, 1):
        if "LHmatrixBlock<5> imnmnmix;" in line:
            linenum = num
            break

with open(temp_location, 'w') as f:
    # Write the stuff at the beginning...
    for i in range(linenum):
        f.write(lines[i])
    # Write the stuff specific to the model.
    f.write((
        "  // LH blocks added by GUM\n"
        "  LHblock<double> dmint;\n"
        "\n"
    ))
    # Then write the rest. Voila: the cheap man's patch.
    for i in range(len(lines)-linenum):
        f.write(lines[i+linenum])

os.remove(location)
os.rename(temp_location, location)
