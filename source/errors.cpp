//    This program, BayeScEnv, aims at detecting genetics markers under local adaptation,
//	based on allele frequency differences between population and environmental differentiation
//    Copyright (C) 2010  Matthieu Foll for the shared code with BayeScan 2.1
//    Copyright (C) 2014  Pierre de Villemereuil for the new code
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.

//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <fstream>
#include <stdlib.h>
//////////////////////////////////
// check file exist
/////////////////////////////////
inline void assure(std::ifstream& in,
                   const char* filename = "")
{
    using namespace std;
    if (!in)
    {
        fprintf(stderr, "Could not open file %s\n", filename);
        exit(EXIT_FAILURE);
    }
}
inline void assure(std::ofstream& in,
                   const char* filename = "")
{
    using namespace std;
    if (!in)
    {
        fprintf(stderr, "Could not open file %s\n", filename);
        exit(EXIT_FAILURE);
    }
}
inline void assure(std::fstream& in,
                   const char* filename = "")
{
    using namespace std;
    if (!in)
    {
        fprintf(stderr, "Could not open file %s\n", filename);
        exit(EXIT_FAILURE);
    }
}



