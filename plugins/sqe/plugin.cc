/*
 * @BEGIN LICENSE
 *
 * sqe by Psi4 Developer, a plugin to:
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include "psi4/psi4-dec.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libpsio/psio.hpp"

#include "sqe.h"

namespace psi{ namespace sqe {

extern "C"
int read_options(std::string name, Options& options)
{
    if (name == "SQE"|| options.read_globals()) {

        /*- the name of the LHS tensor -*/
        options.add_str("SQLHS", "A");
        /*- a string of indices indicating which element of the LHS we are manipulating -*/
        options.add("SQINDICES", new ArrayType());

        /*- a string of creation/annihilation operators to rearrange -*/
        options.add("SQSTRING",new ArrayType());
        /*- a string of indices representing a tensor -*/
        options.add("SQTENSOR",new ArrayType());
        /*- the multiplicative factor for the given string -*/
        options.add_double("SQFACTOR",1.0);

        /*- a string of creation/annihilation operators to rearrange -*/
        options.add("SQSTRING2",new ArrayType());
        /*- a string of indices representing a tensor -*/
        options.add("SQTENSOR2",new ArrayType());
        /*- the multiplicative factor for the given string -*/
        options.add_double("SQFACTOR2",1.0);

        /*- a string of creation/annihilation operators to rearrange -*/
        options.add("SQSTRING3",new ArrayType());
        /*- a string of indices representing a tensor -*/
        options.add("SQTENSOR3",new ArrayType());
        /*- the multiplicative factor for the given string -*/
        options.add_double("SQFACTOR3",1.0);

        /*- a string of creation/annihilation operators to rearrange -*/
        options.add("SQSTRING4",new ArrayType());
        /*- a string of indices representing a tensor -*/
        options.add("SQTENSOR4",new ArrayType());
        /*- the multiplicative factor for the given string -*/
        options.add_double("SQFACTOR4",1.0);
    }

    return true;
}

void removeStar(std::string &x)
{
  auto it = std::remove_if(std::begin(x),std::end(x),[](char c){return (c == '*');});
  x.erase(it, std::end(x));
}

void AddNewString(Options& options,std::vector< ahat* > &ordered,std::string stringnum){

    std::shared_ptr<ahat> mystring (new ahat());

    if ( options["SQFACTOR"+stringnum].has_changed() ) {
        if ( options.get_double("SQFACTOR"+stringnum) > 0.0 ) {
            mystring->sign = 1;
            mystring->factor = fabs(options.get_double("SQFACTOR"+stringnum));
        }else {
            mystring->sign = -1;
            mystring->factor = fabs(options.get_double("SQFACTOR"+stringnum));
        }
    }

    if ( options["SQSTRING"+stringnum].has_changed() ) {
        for (int i = 0; i < (int)options["SQSTRING"+stringnum].size(); i++) {
            std::string me = options["SQSTRING"+stringnum][i].to_string();
            if ( me.find("*") != std::string::npos ) {
                removeStar(me);
                mystring->is_dagger.push_back(true);
            }else {
                mystring->is_dagger.push_back(false);
            }
            mystring->symbol.push_back(me);
        }
    }

    if ( options["SQTENSOR"+stringnum].has_changed() ) {
        for (int i = 0; i < (int)options["SQTENSOR"+stringnum].size(); i++) {
            std::string me = options["SQTENSOR"+stringnum][i].to_string();
            mystring->tensor.push_back(me);
        }
    }

    if ( options["SQLHS"].has_changed() ) {
        mystring->lhs = options["SQLHS"].to_string();
    }

    printf("\n");
    printf("    ");
    printf("// starting string:\n");
    mystring->print();

    // rearrange strings
    mystring->normal_order(ordered);

    // alphabetize
    mystring->alphabetize(ordered);

    // cancel terms
    mystring->cleanup(ordered);

}

void consolidate(std::vector<ahat *> &in,std::vector<ahat *> &out) {

    bool *vanish = (bool*)malloc(in.size()*sizeof(bool));
    memset((void*)vanish,'\0',in.size()*sizeof(bool));
    for (int i = 0; i < (int)in.size(); i++) {
        for (int j = i+1; j < (int)in.size(); j++) {


            bool strings_differ = false;

            // check strings
            if ( in[i]->symbol.size() == in[j]->symbol.size() ) {
                for (int k = 0; k < (int)in[i]->symbol.size(); k++) {

                    // strings differ?
                    if ( in[i]->symbol[k] != in[j]->symbol[k] ) {
                        strings_differ = true;
                    }

                }
            }else {
                strings_differ = true;
            }
            if ( strings_differ ) continue;

            // check deltas
            if ( in[i]->delta1.size() == in[j]->delta1.size() ) {
                for (int k = 0; k < (int)in[i]->delta1.size(); k++) {

                    // strings differ?
                    if ( in[i]->delta1[k] != in[j]->delta1[k] || in[i]->delta2[k] != in[j]->delta2[k] ) {
                        strings_differ = true;
                    }

                }
            }else {
                strings_differ = true;
            }
            if ( strings_differ ) continue;

            // check tensors
            if ( in[i]->tensor.size() == in[j]->tensor.size() ) {
                for (int k = 0; k < (int)in[i]->tensor.size(); k++) {

                    // strings differ?
                    if ( in[i]->tensor[k] != in[j]->tensor[k] ) {

                        strings_differ = true;
                    }

                }
            }else {
                strings_differ = true;
            }
            if ( strings_differ ) continue;
            
            // at this point, we know the strings are the same.  what about the factor?
            int fac1 = in[i]->factor;
            int fac2 = in[j]->factor;
            if ( fabs(fac1 + fac2) < 1e-8 ) {
                vanish[i] = true;
                vanish[j] = true;
                //printf("these terms will cancel\n");
                //in[i]->print();
                //in[j]->print();
            }


        }
    }
    for (int i = 0; i < (int)in.size(); i++) {
        if ( !vanish[i] ) {
            out.push_back(in[i]);
        }
        
    }


}



extern "C"
std::shared_ptr<Wavefunction> sqe(std::shared_ptr<Wavefunction> ref_wfn, Options& options)
{
    std::vector< ahat* > ordered;

    if ( options["SQSTRING"].has_changed() ) {
        AddNewString(options,ordered,"");
    }
    if ( options["SQSTRING2"].has_changed() ) {
        AddNewString(options,ordered,"2");
    }
    if ( options["SQSTRING3"].has_changed() ) {
        AddNewString(options,ordered,"3");
    }
    if ( options["SQSTRING4"].has_changed() ) {
        AddNewString(options,ordered,"4");
    }

    printf("\n");
    printf("    ");
    printf("// normal-ordered strings:\n");
    //for (int i = 0; i < (int)ordered.size(); i++) {
    //    ordered[i]->print();
    //}
    //printf("\n");

    std::vector< ahat* > pruned;
    consolidate(ordered,pruned);

    //printf("\n");
    //printf("some terms cancel:\n");
    for (int i = 0; i < (int)pruned.size(); i++) {
        pruned[i]->check_spin();
        pruned[i]->print();
    }
    printf("\n");

    //printf("%5i %5i\n",(int)ordered.size(),(int)pruned.size());

    // print code?
    for (int i = 0; i < (int)pruned.size(); i++) {
        //pruned[i]->check_spin();
        pruned[i]->print_code(options);
        //pruned[i]->print();
    }
    /*for (int i = 0; i < (int)ordered.size(); i++) {

        std::vector<std::string> sums;
        // check for summation indices
        for (int j = 0; j < (int)ordered[i]->symbol.size(); j++) {
            for (int k = 0; k < (int)ordered[i]->tensor.size(); k++) {
                if ( ordered[i]->symbol[j] == ordered[i]->tensor[k] ) {
                    sums.push_back(ordered[i]->symbol[j]);
                    break;
                }
            }
        }
        for (int j = 0; j < (int)sums.size(); j++) {
            printf("    for (int %s = 0; %s < nmo_; %s++) {\n",sums[j].c_str(),sums[j].c_str(),sums[j].c_str());
        }
    }*/

    return ref_wfn;
}

ahat::ahat() {
}
ahat::~ahat() {
}

void ahat::check_spin() {

    // check A/B in delta functions
    for (int j = 0; j < (int)delta1.size(); j++) {
        if ( delta1[j].length() == 2 ) {
            if ( delta1[j].at(1) == 'A' && delta2[j].at(1) == 'B' ) {
                skip = true;
                break;
            }else if ( delta1[j].at(1) == 'B' && delta2[j].at(1) == 'A' ) {
                skip = true;
                break;
            }
        }
    }

    // check A/B in two-index tensors
    if ( (int)tensor.size() == 2 ) {
        if ( tensor[0].length() == 2 ) {
            if ( tensor[1].length() == 2 ) {

                if ( tensor[0].at(1) == 'A' && tensor[1].at(1) == 'B' ) {
                    skip = true;
                    return;
                }else if ( tensor[0].at(1) == 'B' && tensor[1].at(1) == 'A' ) {
                    skip = true;
                    return;
                }
            
            }
        }
    }

    // check A/B in four-index tensors
    if ( (int)tensor.size() == 4 ) {
        // check bra
        if ( tensor[0].length() == 2 ) {
            if ( tensor[1].length() == 2 ) {

                if ( tensor[0].at(1) == 'A' && tensor[1].at(1) == 'B' ) {
                    skip = true;
                    return;
                }else if ( tensor[0].at(1) == 'B' && tensor[1].at(1) == 'A' ) {
                    skip = true;
                    return;
                }
            
            }
        }
        // check ket
        if ( tensor[2].length() == 2 ) {
            if ( tensor[3].length() == 2 ) {

                if ( tensor[2].at(1) == 'A' && tensor[3].at(1) == 'B' ) {
                    skip = true;
                    return;
                }else if ( tensor[2].at(1) == 'B' && tensor[3].at(1) == 'A' ) {
                    skip = true;
                    return;
                }
            
            }
        }

    }


}

void ahat::print_code(Options& options) {
    if ( skip ) return;

    std::vector<std::string> sums;
    // check for summation indices
    for (int j = 0; j < (int)symbol.size(); j++) {
        for (int k = 0; k < (int)tensor.size(); k++) {
            if ( symbol[j] == tensor[k] ) {
                sums.push_back(symbol[j]);
                break;
            }
        }
    }

    int indent = 1;

    // first set of indices corresponds to normal-ordered operator
    for (int j = 0; j < (int)symbol.size(); j++) {
        bool skipme = false;
        for (int k = 0; k < (int)sums.size(); k++) {

            // skip this index if it is a summation index
            if ( symbol[j] == sums[k] ) {
                skipme = true;
                break;
            }
        }
        if ( skipme ) continue;
        for (int k = 0; k < indent; k++) {
            printf("    ");
        }
        printf("for (int %s = 0; %s < nmo_; %s++) {\n",symbol[j].c_str(),symbol[j].c_str(),symbol[j].c_str());

        indent++;
    }
    // second set of indices corresponds to other tensor
    for (int j = 0; j < (int)tensor.size(); j++) {
        bool skipme = false;
        for (int k = 0; k < (int)sums.size(); k++) {

            // skipme this index if it is a summation index
            if ( tensor[j] == sums[k] ) {
                skipme = true;
                break;
            }
        }
        if ( skipme ) continue;
        for (int k = 0; k < indent; k++) {
            printf("    ");
        }
        printf("for (int %s = 0; %s < nmo_; %s++) {\n",tensor[j].c_str(),tensor[j].c_str(),tensor[j].c_str());

        indent++;
    }

    // third set of indices accounts for any leftover indices on lhs that arise in delta functions only

    std::vector<std::string> delta_loops;
    for (int j = 0; j < (int)delta1.size(); j++) {
        for (int k = 0; k < indent; k++) {
            printf("    ");
        }
        printf("for (int %s = 0; %s < nmo_; %s++) {\n",delta1[j].c_str(),delta1[j].c_str(),delta1[j].c_str());

        indent++;
    }


    for (int k = 0; k < indent; k++) {
        printf("    ");
    }
    printf("double dum = 0.0;\n");
    for (int j = 0; j < (int)sums.size(); j++) {
        for (int k = 0; k < indent; k++) {
            printf("    ");
        }
        printf("for (int %s = 0; %s < nmo_; %s++) {\n",sums[j].c_str(),sums[j].c_str(),sums[j].c_str());

        indent++;
    }

    // LHS[INDICES] += FACTOR * SIGN * TENSOR * OPERATOR

    //printf("%s[] += %20.12le * \n",lhs.c_str(),factor*sign);
    std::string tensor_string;
    if ( (int)tensor.size() == 4 ) {
        tensor_string.append("tei[INDEX(");
        tensor_string.append(tensor[0].c_str());
        tensor_string.append(",");
        tensor_string.append(tensor[1].c_str());
        tensor_string.append(")*nmo_*(nmo+1)/2+INDEX(");
        tensor_string.append(tensor[2].c_str());
        tensor_string.append(",");
        tensor_string.append(tensor[3].c_str());
        tensor_string.append(")]");
    }else if ( (int)tensor.size() == 2 ) {
        tensor_string.append("h1[");
        tensor_string.append(tensor[0].c_str());
        tensor_string.append("*nmo_+");
        tensor_string.append(tensor[1].c_str());
        tensor_string.append("]");
    }

    std::string operator_string;
    if ( symbol.size() == 4 ) {

        int na = 0;
        int nb = 0;
        if ( symbol[0].size() != 2 || symbol[1].size() != 2 || symbol[2].size() != 2 || symbol[3].size() != 2 ) {
            throw PsiException("SQE only works with spin labels",__FILE__,__LINE__);
        }
        if ( symbol[0].at(1) == 'A' ) na++;
        else                          nb++;
        if ( symbol[1].at(1) == 'A' ) na++;
        else                          nb++;
        if ( symbol[2].at(1) == 'A' ) na++;
        else                          nb++;
        if ( symbol[3].at(1) == 'A' ) na++;
        else                          nb++;

        if ( na % 2 != 0 || nb % 2 != 0 ) {
            throw PsiException("something is wrong. spin symmetry broken in ordered operator.",__FILE__,__LINE__);
        }

        if ( na == 4 ) {
            operator_string.append("D2aa[");
            operator_string.append(symbol[0]);
            operator_string.append("*nmo_*nmo_*nmo_+");
            operator_string.append(symbol[1]);
            operator_string.append("*nmo_*nmo_+");
            operator_string.append(symbol[3]);
            operator_string.append("*nmo_+");
            operator_string.append(symbol[2]);
            operator_string.append("]");
        }else if ( nb == 4 ) {
            operator_string.append("D2bb[");
            operator_string.append(symbol[0]);
            operator_string.append("*nmo_*nmo_*nmo_+");
            operator_string.append(symbol[1]);
            operator_string.append("*nmo_*nmo_+");
            operator_string.append(symbol[3]);
            operator_string.append("*nmo_+");
            operator_string.append(symbol[2]);
            operator_string.append("]");
        } else {
            if ( symbol[0].at(1) == 'B' && symbol[1].at(1) == 'A' ) {
                sign = -sign;
                if ( symbol[3].at(1) == 'B' && symbol[2].at(1) == 'A' ) {
                    sign = -sign;
                    operator_string.append("D2ab[");
                    operator_string.append(symbol[1]);
                    operator_string.append("*nmo_*nmo_*nmo_+");
                    operator_string.append(symbol[0]);
                    operator_string.append("*nmo_*nmo_+");
                    operator_string.append(symbol[2]);
                    operator_string.append("*nmo_+");
                    operator_string.append(symbol[3]);
                    operator_string.append("]");
                }else {
                    operator_string.append("D2ab[");
                    operator_string.append(symbol[1]);
                    operator_string.append("*nmo_*nmo_*nmo_+");
                    operator_string.append(symbol[0]);
                    operator_string.append("*nmo_*nmo_+");
                    operator_string.append(symbol[3]);
                    operator_string.append("*nmo_+");
                    operator_string.append(symbol[2]);
                    operator_string.append("]");
                }
            }else {
                if ( symbol[3].at(1) == 'B' && symbol[2].at(1) == 'A' ) {
                    sign = -sign;
                    operator_string.append("D2ab[");
                    operator_string.append(symbol[0]);
                    operator_string.append("*nmo_*nmo_*nmo_+");
                    operator_string.append(symbol[1]);
                    operator_string.append("*nmo_*nmo_+");
                    operator_string.append(symbol[2]);
                    operator_string.append("*nmo_+");
                    operator_string.append(symbol[3]);
                    operator_string.append("]");
                }else {
                    operator_string.append("D2ab[");
                    operator_string.append(symbol[0]);
                    operator_string.append("*nmo_*nmo_*nmo_+");
                    operator_string.append(symbol[1]);
                    operator_string.append("*nmo_*nmo_+");
                    operator_string.append(symbol[3]);
                    operator_string.append("*nmo_+");
                    operator_string.append(symbol[2]);
                    operator_string.append("]");
                }
            }
        }

    }else if ( symbol.size() == 2 ) {

        int na = 0;
        int nb = 0;
        if ( symbol[0].size() != 2 || symbol[1].size() != 2 ) {
            throw PsiException("SQE only works with spin labels",__FILE__,__LINE__);
        }
        if ( symbol[0].at(1) == 'A' ) na++;
        else                          nb++;
        if ( symbol[1].at(1) == 'A' ) na++;
        else                          nb++;

        if ( na % 2 != 0 || nb % 2 != 0 ) {
            throw PsiException("something is wrong. spin symmetry broken in ordered operator.",__FILE__,__LINE__);
        }

        if ( na == 2 ) {

            operator_string.append("D1a[");
            operator_string.append(symbol[0]);
            operator_string.append("*nmo_+");
            operator_string.append(symbol[1]);
            operator_string.append("]");

        }else {

            operator_string.append("D1b[");
            operator_string.append(symbol[0]);
            operator_string.append("*nmo_+");
            operator_string.append(symbol[1]);
            operator_string.append("]");

        }


    }else {
        printf("// warning, skipping contributions from %i-particle RDMs.\n",symbol.size()/2);
    }

    if ( symbol.size() == 4 || symbol.size() == 2 ) {

        for (int k = 0; k < indent; k++) {
            printf("    ");
        }
        printf("dum += %20.12le\n",factor*sign);
            
        if ( (int)tensor.size() > 0 ) {
            for (int k = 0; k < indent; k++) {
                printf("    ");
            }
            printf("     * %s\n",tensor_string.c_str());
        }

        for (int k = 0; k < indent; k++) {
            printf("    ");
        }
        printf("     * %s;\n",operator_string.c_str());
            
    }

    // one too many indents
    indent--;

    // close brackets for summation indices
    for (int i = 0; i < sums.size(); i++) {
        for (int j = 0; j < indent; j++) {
            printf("    ");
        }
        printf("}\n");
        indent--;
    }
    printf("\n");

    // accumulate result into lhs
    std::vector<std::string> indices;
    for (int i = 0; i < (int)options["SQINDICES"].size(); i++) {
        std::string me = options["SQINDICES"][i].to_string();
        indices.push_back(me);
    }

    // account for additional delta functions
    for (int i = 0; i < (int)delta1.size(); i++) {
        for (int j = 0; j < (int)indices.size(); j++) {
            if ( delta2[i] == indices[j] ) {
                indices[j] = delta1[i];
            }
        }
    }

    for (int j = 0; j < indent; j++) {
        printf("    ");
    }
    printf("    ");

    if ( indices.size() == 4 ) {
    
        if ( indices[0].at(1) != indices[1].at(1) ) {
            throw PsiException("sqe only supports spin-conserving tensors at this point",__FILE__,__LINE__);
        }else if ( indices[2].at(1) != indices[3].at(1) ) {
            throw PsiException("sqe only supports spin-conserving tensors at this point",__FILE__,__LINE__);
        }

        // if aa: (%s*nmo_+%s            )
        // if bb: (%s*nmo_+%s + nmo_*nmo_)
        std::string bra = "(";
        bra.append(indices[0]);
        bra.append("*nmo_+");
        bra.append(indices[1]);
        if ( indices[0].at(1) == 'A' ) {
            bra.append("            )");
        }else {
            bra.append(" + nmo_*nmo_)");
        }

        std::string ket = "(";
        ket.append(indices[2]);
        ket.append("*nmo_+");
        ket.append(indices[3]);
        if ( indices[2].at(1) == 'A' ) {
            ket.append("            )");
        }else {
            ket.append(" + nmo_*nmo_)");
        }


        printf("%s[%s*2*nmo_*nmo_+%s] += dum;\n",
            lhs.c_str(),
            bra.c_str(),
            ket.c_str());
        //printf("%s[%s*nmo_*nmo_*nmo_+%s*nmo_*nmo_+%s*nmo_+%s] += dum;\n",
        //    lhs.c_str(),
        //    indices[0].c_str(),
        //    indices[1].c_str(),
        //    indices[2].c_str(),
        //    indices[3].c_str());
    } else if ( indices.size() == 2 ) {
        //printf("%s[%s*nmo_+%s] += dum;\n",
        //    lhs.c_str(),
        //    indices[0].c_str(),
        //    indices[1].c_str());
        throw PsiException("i don't understand the dimensions of the tensor on the lhs.",__FILE__,__LINE__);
    } else {
        throw PsiException("i don't understand the dimensions of the tensor on the lhs.",__FILE__,__LINE__);
    }

    // close brackets for remaining indices
    for (int i = 0; i < indices.size() - delta1.size(); i++) {
        for (int j = 0; j < indent; j++) {
            printf("    ");
        }
        printf("}\n");
        indent--;
    }
    printf("\n");
}

void ahat::print() {
    if ( skip ) return;
    printf("    ");
    printf("//     ");
    printf("%c", sign > 0 ? '+' : '-');
    printf(" ");
    printf("%4.2lf", fabs(factor));
    printf(" ");
    for (int i = 0; i < (int)symbol.size(); i++) {
        printf("%s",symbol[i].c_str());
        if ( is_dagger[i] ) {
            printf("%c",'*');
        }
        printf(" ");
    }
    for (int i = 0; i < (int)delta1.size(); i++) {
        printf("d(%s%s)",delta1[i].c_str(),delta2[i].c_str());
        printf(" ");
    }
    if ( (int)tensor.size() > 0 ) {
        // two-electron integrals
        if ( (int)tensor.size() == 4 ) {
            printf("(");
            for (int i = 0; i < 2; i++) {
                printf("%s",tensor[i].c_str());
            }
            printf("|");
            for (int i = 2; i < 4; i++) {
                printf("%s",tensor[i].c_str());
            }
            printf(")");
        }
        // one-electron integrals
        if ( (int)tensor.size() == 2 ) {
            printf("h(");
            for (int i = 0; i < 2; i++) {
                printf("%s",tensor[i].c_str());
            }
            printf(")");
        }
    }
    printf("\n");
}

bool ahat::is_normal_order() {

    // don't bother bringing to normal order if we're going to skip this string
    if (skip) return true;

    for (int i = 0; i < (int)symbol.size()-1; i++) {
        if ( !is_dagger[i] && is_dagger[i+1] ) {
            return false;
        }
    }
    return true;
}

// in order to compare strings, the creation and annihilation 
// operators should be ordered in some consistent way.
// alphabetically seems reasonable enough
void ahat::alphabetize(std::vector<ahat *> &ordered) {

    // alphabetize string
    for (int i = 0; i < (int)ordered.size(); i++) {

        // creation
        bool not_alphabetized = false;
        do {
            not_alphabetized = false;
            int ndagger = 0;
            for (int j = 0; j < (int)ordered[i]->symbol.size(); j++) {
                if ( ordered[i]->is_dagger[j] ) ndagger++;
            }
            for (int j = 0; j < ndagger-1; j++) {
                int val1 = ordered[i]->symbol[j].c_str()[0];
                int val2 = ordered[i]->symbol[j+1].c_str()[0];
                if ( val2 < val1 ) {
                    std::string dum = ordered[i]->symbol[j];
                    ordered[i]->symbol[j] = ordered[i]->symbol[j+1];
                    ordered[i]->symbol[j+1] = dum;
                    ordered[i]->sign = -ordered[i]->sign;
                    not_alphabetized = true;
                    j = (int)ordered[i]->symbol.size() + 1;
                    not_alphabetized = true;
                }
            }
        }while(not_alphabetized);
        // annihilation
        not_alphabetized = false;
        do {
            not_alphabetized = false;
            int ndagger = 0;
            for (int j = 0; j < (int)ordered[i]->symbol.size(); j++) {
                if ( ordered[i]->is_dagger[j] ) ndagger++;
            }
            for (int j = ndagger; j < (int)ordered[i]->symbol.size()-1; j++) {
                int val1 = ordered[i]->symbol[j].c_str()[0];
                int val2 = ordered[i]->symbol[j+1].c_str()[0];
                if ( val2 < val1 ) {
                    std::string dum = ordered[i]->symbol[j];
                    ordered[i]->symbol[j] = ordered[i]->symbol[j+1];
                    ordered[i]->symbol[j+1] = dum;
                    ordered[i]->sign = -ordered[i]->sign;
                    not_alphabetized = true;
                    j = (int)ordered[i]->symbol.size() + 1;
                    not_alphabetized = true;
                }
            }
        }while(not_alphabetized);
    }

    // alphabetize deltas
    for (int i = 0; i < (int)ordered.size(); i++) {
        for (int j = 0; j < (int)ordered[i]->delta1.size(); j++) {
            int val1 = ordered[i]->delta1[j].c_str()[0];
            int val2 = ordered[i]->delta2[j].c_str()[0];
            if ( val2 < val1 ) {
                std::string dum = ordered[i]->delta1[j];
                ordered[i]->delta1[j] = ordered[i]->delta2[j];
                ordered[i]->delta2[j] = dum;
            }
        }
    }
}

// once strings are alphabetized, we can compare them
// and remove terms that cancel
void ahat::cleanup(std::vector<ahat *> &ordered) {

    for (int i = 0; i < (int)ordered.size(); i++) {

        for (int j = i+1; j < (int)ordered.size(); j++) {
            
            // same factor
            if ( ordered[i]->factor == ordered[j]->factor ) {

                // opposite sign
                if ( ordered[i]->sign == -ordered[j]->sign ) {

                    // same normal-ordered operator
                    if ( ordered[i]->symbol.size() == ordered[j]->symbol.size() ) {
                        int nsame_s = 0;
                        for (int k = 0; k < (int)ordered[i]->symbol.size(); k++) {
                            if ( ordered[i]->symbol[k] == ordered[j]->symbol[k] ) {
                                nsame_s++;
                            }
                        }
                        if ( nsame_s == ordered[i]->symbol.size() ) {
                            // same tensor
                            if ( ordered[i]->tensor.size() == ordered[j]->tensor.size() ) {
                                int nsame_t = 0;
                                for (int k = 0; k < (int)tensor.size(); k++) {
                                    if ( ordered[i]->tensor[k] == ordered[j]->tensor[k] ) {
                                        nsame_t++;
                                    }
                                }
                                if ( nsame_t == ordered[i]->tensor.size() ) {
                                    // same delta functions (recall these aren't sorted in any way)
                                    int nsame_d = 0;
                                    for (int k = 0; k < (int)ordered[i]->delta1.size(); k++) {
                                        for (int l = 0; l < (int)ordered[j]->delta1.size(); l++) {
                                            if ( ordered[i]->delta1[k] == ordered[j]->delta1[l] && ordered[i]->delta2[k] == ordered[j]->delta2[l] ) {
                                                nsame_d++;
                                            }else if ( ordered[i]->delta2[k] == ordered[j]->delta1[l] && ordered[i]->delta1[k] == ordered[j]->delta2[l] ) {
                                                nsame_d++;
                                            }
                                        }
                                    }
                                    if ( nsame_d == (int)ordered[i]->delta1.size() ) {
                                        ordered[i]->skip = true;
                                        ordered[j]->skip = true;
                                    }
                                }
                            }
                        }
                    }

                }

            }
            
        }

    }

}

void ahat::normal_order(std::vector<ahat *> &ordered) {
    if ( skip ) return;

    if ( is_normal_order() ) {

        // push current ordered operator onto running list
        ahat * newguy (new ahat());

        newguy->skip   = skip;
        newguy->lhs    = lhs;
        newguy->sign   = sign;
        newguy->factor = factor;
        for (int j = 0; j < (int)is_dagger.size(); j++) {
            newguy->is_dagger.push_back(is_dagger[j]);
        }
        for (int j = 0; j < (int)symbol.size(); j++) {
            newguy->symbol.push_back(symbol[j]);
        }
        for (int j = 0; j < (int)tensor.size(); j++) {
            // does tensor index show up in a delta function?
            bool skipme = false;
            for (int k = 0; k < (int)delta1.size(); k++) {
                if ( tensor[j] == delta1[k] ) {
                    newguy->tensor.push_back(delta2[k]);
                    skipme = true;
                    break;
                }
                if ( tensor[j] == delta2[k] ) {
                    newguy->tensor.push_back(delta1[k]);
                    skipme = true;
                    break;
                }
            }
            if ( skipme ) continue;
            newguy->tensor.push_back(tensor[j]);
        }
        for (int j = 0; j < (int)delta1.size(); j++) {
            bool skipme = false;
            for (int k = 0; k < (int)tensor.size(); k++) {
                if ( tensor[k] == delta1[j] ) {
                    skipme = true;
                    break;
                }
                if ( tensor[k] == delta2[j] ) {
                    skipme = true;
                    break;
                }
            }
            if ( skipme ) continue;

            newguy->delta1.push_back(delta1[j]);
            newguy->delta2.push_back(delta2[j]);
        }

        ordered.push_back(newguy);

        return;
    }

    // new strings
    std::shared_ptr<ahat> s1 ( new ahat() );
    std::shared_ptr<ahat> s2 ( new ahat() );

    for (int i = 0; i < (int)tensor.size(); i++) {
        s1->tensor.push_back(tensor[i]);
        s2->tensor.push_back(tensor[i]);
    }

    s1->skip = skip;
    s2->skip = skip;

    s1->sign = sign;
    s2->sign = sign;

    s1->factor = factor;
    s2->factor = factor;

    s1->lhs = lhs;
    s2->lhs = lhs;

    for (int i = 0; i < (int)delta1.size(); i++) {
        s1->delta1.push_back(delta1[i]);
        s2->delta1.push_back(delta1[i]);

        s1->delta2.push_back(delta2[i]);
        s2->delta2.push_back(delta2[i]);
    }

    for (int i = 0; i < (int)symbol.size()-1; i++) {

        if ( !is_dagger[i] && is_dagger[i+1] ) {

            s1->delta1.push_back(symbol[i]);
            s1->delta2.push_back(symbol[i+1]);

            // check spin in delta functions
            for (int j = 0; j < (int)delta1.size(); j++) {
                if ( s1->delta1[j].length() != 2 ) {
                    throw PsiException("be sure to specify spin as second character in labels",__FILE__,__LINE__);
                }
                if ( s1->delta1[j].at(1) == 'A' && s1->delta2[j].at(1) == 'B' ) {
                    s1->skip = true;
                }else if ( s1->delta1[j].at(1) == 'B' && s1->delta2[j].at(1) == 'A' ) {
                    s1->skip = true;
                }
            }

            s2->sign = -s2->sign;
            s2->symbol.push_back(symbol[i+1]);
            s2->symbol.push_back(symbol[i]);
            s2->is_dagger.push_back(is_dagger[i+1]);
            s2->is_dagger.push_back(is_dagger[i]);

            for (int j = i+2; j < (int)symbol.size(); j++) {

                s1->symbol.push_back(symbol[j]);
                s2->symbol.push_back(symbol[j]);

                s1->is_dagger.push_back(is_dagger[j]);
                s2->is_dagger.push_back(is_dagger[j]);

            }
            break;

        }else {

            s1->symbol.push_back(symbol[i]);
            s2->symbol.push_back(symbol[i]);

            s1->is_dagger.push_back(is_dagger[i]);
            s2->is_dagger.push_back(is_dagger[i]);

        }
    }

    s1->normal_order(ordered);
    s2->normal_order(ordered);

}

}} // End namespaces

