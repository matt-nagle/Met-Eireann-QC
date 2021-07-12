#include <vector>
#include <math.h>
#include "titanlibcustom.h"
#include <assert.h>
#include <iostream>
#include <numeric>
#include <exception>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_linalg.h>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>

using namespace titanlibcustom;

//=============================================================================

// FGT within the inner circle

bool fgt_core( const vec& lats, const vec& lons, const vec& yo, const vec& yb, const vec& sigma_b, float minp, float maxp, const vec& mina, const vec& maxa, const vec& minv, const vec& maxv, const vec& tpos, const vec& tneg, const ivec& indices_global_test, const ivec& indices_outer_inner, const ivec& indices_outer_test, const ivec& indices_inner_test, bool debug, float na, bool set_flag0, int& thrown_out, vec& scores, ivec& flags);


//=============================================================================

//+ First guess test 
ivec titanlibcustom::fgt(const Points& points,
                    const vec& values,
                    const ivec& obs_to_check,
                    const vec& background_values,
                    const vec& background_uncertainties,
                    BackgroundType background_elab_type,
                    int num_min_outer,
                    int num_max_outer,
                    float inner_radius,
                    float outer_radius,
                    int num_iterations,
                    int num_min_prof,
                    float min_elev_diff,
                    const vec& values_mina,
                    const vec& values_maxa,
                    const vec& values_minv,
                    const vec& values_maxv,
                    const vec& tpos,
                    const vec& tneg,
                    bool debug,
                    vec& scores) {
/*

 -+-+ First guess test +-+-

 + Description: 
    Flag observations that are (likely) affected by gross measurement errors (GE) based on neighbouring observations. Observations are affected by GEs if their values are (a) [not related to the actual atmospheric state] OR (b) [affected by such large representativeness errors (REs) that they are difficult to reconstruct using the neighbouring observations].
  For simplicity, we will refer to observations affected by GE as bad observations, viceversa good observations are those not affected by GEs. 
  FGT is a simplified SCT.

 + Algorithm: 
    The big task of applying the FGT over the whole observational network is splitted into several smaller tasks.
  The FGT is applied many times in sequence. 
  Each time, the FGT is applied over a small region centered on an observation location (i.e. the centroid observation).
  Based on the centroid, we define an inner circle and an outer circle.
  The outer circle is used to obtain a spatial trend of the observed variable over the region, the so-called background.
  There are many ways to compute such background values.
  In addition, the backgreund can be provided as an input vector, together with its uncertainty.
  The observations to test are those within the inner circle that have not been tested yet.
  See also the comments of the function "sct_core" for more details.
  
  The application of the sequence of FGTs over the whole observational network is repeated several times, until no GEs are found within the observations tested. 
  The core of FGT is to test simultaneously some observations assuming (almost) all the neighbouring observations are good ones. As a consequence, bad observations can influence the results of some tests on neighbouring observations before being rejected.
  We apply special care to avoid "false positives" (good flagged as bad).
  For each FGT over a centroid, we flag as bad only the observation showing the worst result.
  In addition, an extra-loop is done over bad observations where only the good ones are considered. We may "save" bad observations.
  Analogously, to avoid "false negatives" (bad flagged as good), the last loop is done over the bad observations where only the good ones are considered. We may "save" bad observations during this last loop.

  (*) Notes on prior information (see also comments in sct_core): the user is expected to provide two different range of values, that inform about the desired tolerance level of the FGT.
    i.  The range of admissible values (values_mina,values_maxa). One range for each observation. The observed value have to be within the corresponding range (e.g. observation=2, then OK if values_mina=0 & values_maxa=4; not OK if values_mina=10 & values_maxa=14). A background outside this range warns us about possible problems with the observation.
        a.   The range of plausible values (value_minp,value_maxp) is derived from the range of admissible values and it is used to check the analyses and backgrounds. Outside this range the values do not make sense.
    ii. The range of valid values (values_minv,values_maxv). One range for each observation. The observed value must be inside this range. If the backgrounds in a region are all within this range, then we assume the corresponding observations are good ones.
  Examples:
    Temperature: 3 observed values(degC)=( 10, 15, 5). (values_mina,values_maxa) = observed values plus/minus 15 degC. (values_minv,values_maxv) = observed values plus/minus 1 degC.
    Precipitation: 3 observed values(mm)=( 0, 50, 1). (values_mina,values_maxa) = observed values plus/minus either 10 mm or 50% of the value. (values_minv,values_maxv) = observed values plus/minus either 1 mm or 10% of the value.

  (*) Notes on the data structures (see also comments in set_indices): we have defined a hierarchy of data structures
    global -> outer -> inner -> test
   global vectors are the input vectors. The other data structures refer to a specific iteration when a centroid observation has been identified. outer are the vectors of data on the outer circle, not yet flagged as bad. inner are the data on the inner circle, not yet flagged as bad. test data are the observations in the inner circle and still to test.
   We use the same counter variables for the data structures: g -> global; i,j -> outer; l -> inner; m -> test
   There are vectors of indices matching the different structures. 
   For example, indices_outer_inner matches the outer and the inner vectors. indices_outer_inner is a vector of the dimension of the number of observations in the inner circle. The instruction i=indices_outer_inner[l] maps the l-th element of inner-type vectors as the i-th element of outer-type vectors.

 (*) Notes on the thresholds. Different thresholds can be used for different observations. The thresholds vectors have the same dimensions of the observations vectors (global). For each observation, two thresholds (tpos,tneg) have to be set by the user. tpos is the threshold used when the observed value is higher than the cv-analysis, otherwise tneg is used. 

 + Pseudo-code:
 
  Loop over the number of iterations prescribed
    Loop over the p observations (observation index is "curr"). 
      Decide if the observation is suitable for a centroid
      Define outer/inner cirles, observations to test and set all the vector of indices
      Check if the observation is isolated
      Compute the background
      If the deviations between observations and background is small enough that further steps are not needed and the observations to test are flagged as good ones
      FGT over the inner circle, either flag the worse observation as bad or flag all of them as good ones
      NOTE: observations can be flagged as good ones only after the first iteration (after rejecting the worst)
  
  Loop over the observations without QC flags (nor good neither bad) and use each of them as a centroid

  Loop over the bad observations and use each of them as a centroid
  
 + Returned values:

  scores. only for observations rejected by the FGT, z-score (See fgt_core)

  flags. -999 = not checked; 0 = passed (good); 1 = failed (bad); 11 = isolated (<2 inside inner); 12 = isolated (<num_min_outer inside outer)

*/
    double s_time = titanlibcustom::clock();
    const vec& lats = points.get_lats();
    const vec& lons = points.get_lons();
    const vec& elevs = points.get_elevs();

    const int p = values.size();
    if( lats.size() != p || lons.size() != p || elevs.size() != p || values.size() != p || tpos.size() != p || tneg.size() != p || values_mina.size() != p || values_maxa.size() != p || values_minv.size() != p || values_maxv.size() != p)
        throw std::runtime_error("Dimension mismatch");
    if(num_min_outer < 2)
        throw std::invalid_argument("num_min_outer must be > 1");
    if(num_max_outer < num_min_outer)
        throw std::invalid_argument("num_max_outer must be > num_min_outer");
    if(num_iterations < 1)
        throw std::invalid_argument("num_iterations must be >= 1");
    if(min_elev_diff <= 0)
        throw std::invalid_argument("min_elev_diff must be > 0");
    if(inner_radius < 0)
        throw std::invalid_argument("inner_radius must be >= 0");
    if(outer_radius < inner_radius)
        throw std::invalid_argument("outer_radius must be >= inner_radius");
    if(background_elab_type == titanlibcustom::VerticalProfile && num_min_prof<0)
        throw std::invalid_argument("num_min_prof must be >=0");
    if(background_elab_type == titanlibcustom::External &&  background_values.size() != p)
        throw std::runtime_error("Background vector dimension mismatch");

    // initializations
    float na = -999.; // code for "Not Available". Any type of missing data
    ivec flags( p, na);
    ivec obs_test( p, 1);
    vec sigma_b( p, 1);

    bool set_all_good = false;
    
    scores.clear();
    scores.resize( p, na);
    
    if (obs_to_check.size() == p) {
        for(int g=0; g < p; g++) 
          obs_test[g] = obs_to_check[g];
    }
 
    if (background_uncertainties.size() == p) {
        for(int g=0; g < p; g++) 
          sigma_b[g] = background_uncertainties[g];
    }

    /* plausible values are inferred from admissible values
       note: a range check should be performed before doing the FGT.
             The FGT can mistake obvious large errors for good observations 
             if they happen to be close enough that they support each other */
    float value_minp = values_mina[0];
    float value_maxp = values_maxa[0];
    for(int g=1; g < p; g++) {
        if ( values_mina[g] < value_minp) value_minp = values_mina[g];
        if ( values_maxa[g] > value_maxp) value_maxp = values_maxa[g];
    }
    if(debug) std::cout << "=============================================== " << p << std::endl;
    if(debug) std::cout << "Number of observations to test is " << p << std::endl;
    if(debug) std::cout << "Min and max acceptable values are " << value_minp << " " << value_maxp << std::endl;

    // KDtree has to do with fast computation of distances
    titanlibcustom::KDTree tree(lats, lons);

    //
    if(debug) {
        std::cout << "g lats lons obs:" << std::endl;
        for(int g=0; g < p; g++) {
            std::cout << g << " " << lats[g] << " " << lons[g] << " " << values[g] << std::endl;
        }
    }

    // FGT iterations
    for(int iteration = 0; iteration < num_iterations; iteration++) {

        if(debug) std::cout << " +++++ Iteration " << iteration << " ++++++++++++++++" << std::endl;

        double s_time0 = titanlibcustom::clock();
        
        // reset this number each loop (this is for breaking if we don't throw anything new out)
        int thrown_out = 0; 
        
        // diagnostic. count the number of times the test is performed
        int count_tests = 0;

        // loop over observations
        for(int curr=0; curr < p; curr++) {

            if(debug) std::cout << "===> curr " << curr << " ===============" << std::endl;

            // -~- Can the observation be used as a centroid for the definition of inner and outer circles?

            // jump to next if observation not to be checked or already checked 
            if(obs_test[curr] != 1 || flags[curr] >= 0) {
                if(debug) std::cout << "..observation not suitable as a centroid " << curr << std::endl;
                continue;
            }

            // -~- Define outer/inner circles and which are the observations to test

            // get all neighbours that are close enough (inside outer circle)
            vec distances;
            ivec indices_global_outer_guess = tree.get_neighbours_with_distance(lats[curr], lons[curr], outer_radius, distances, true);
            if(indices_global_outer_guess.size() > num_max_outer) {

                int N = indices_global_outer_guess.size();
                std::vector<std::pair<float,int> > pairs(N);
                for(int i = 0; i < indices_global_outer_guess.size(); i++) {
                    pairs[i] = std::pair<float, int>(distances[i], indices_global_outer_guess[i]);
                }
                std::sort(pairs.begin(), pairs.end(), titanlibcustom::sort_pair_first<float,int>());
                distances.clear();
                indices_global_outer_guess.clear();
                distances.resize(num_max_outer);
                indices_global_outer_guess.resize(num_max_outer);
                for(int i = 0; i < num_max_outer; i++) {
                    distances[i] = pairs[i].first;
                    indices_global_outer_guess[i] = pairs[i].second;
                }
            }
            if(debug) {
                int p_dist = distances.size();
                std::cout << "p_dist  " << p_dist << std::endl;
            }

            // set all the indices linking the different levels

            ivec indices_global_outer, indices_global_test, indices_outer_inner, indices_outer_test, indices_inner_test;
            bool res = titanlibcustom::set_indices(indices_global_outer_guess, obs_test, flags, distances, inner_radius, -1, indices_global_outer, indices_global_test, indices_outer_inner, indices_outer_test, indices_inner_test);

            int p_outer = indices_global_outer.size();
            int p_inner = indices_outer_inner.size();
            int p_test = indices_global_test.size();
            if(debug) std::cout << "p_outer inner test " << p_outer << " " << p_inner << " " << "" << p_test << std::endl;

            // -~- Check if there are enough observations in the outer/inner circles

            if(p_outer < num_min_outer) {
                flags[curr] = 12;
                if(debug) std::cout << "@@isolated (outer) " << curr << std::endl;
                continue; 
            }

            if( p_inner < 2) {
                flags[curr] = 11;
                if(debug) std::cout << "@@isolated (inner) " << curr << std::endl;
                continue;
            }
            
            // -~- Vectors on the outer circle

            vec lons_outer   = titanlibcustom::subset(lons, indices_global_outer);
            vec elevs_outer  = titanlibcustom::subset(elevs, indices_global_outer);
            vec lats_outer   = titanlibcustom::subset(lats, indices_global_outer);
            vec sigmab_outer = titanlibcustom::subset(sigma_b, indices_global_outer);
            vec values_outer = titanlibcustom::subset(values, indices_global_outer);
            vec tpos_outer   = titanlibcustom::subset(tpos, indices_global_outer);
            vec tneg_outer   = titanlibcustom::subset(tneg, indices_global_outer);
            vec mina_outer   = titanlibcustom::subset(values_mina, indices_global_outer);
            vec minv_outer   = titanlibcustom::subset(values_minv, indices_global_outer);
            vec maxa_outer   = titanlibcustom::subset(values_maxa, indices_global_outer);
            vec maxv_outer   = titanlibcustom::subset(values_maxv, indices_global_outer);
            
            // Debug for indices
            if(debug) {
                ivec flags_outer;
                flags_outer.resize( p_outer, na);
                for(int i=0; i < p_outer; i++) {
                   int g = indices_global_outer[i];
                   flags_outer[i] = flags[g];
                }
                std::cout << "curr lats lons obs:" << std::endl;
                std::cout << curr << " " << lats[curr] << " " << lons[curr] << " " << values[curr] << std::endl;
                std::cout << "indices_global_outer - i g lats lons obs:" << std::endl;
                for(int i=0; i<p_outer; i++) {
                    int g = indices_global_outer[i];
                    std::cout << std::setprecision(6) << i << " " << g << " " << lats[g] << " " << lons[g] << " " << values[g] << " " << flags[g] << std::endl;
                }
                std::cout << "outer - i lats lons obs:" << std::endl;
                for(int i=0; i<p_outer; i++) {
                    std::cout << std::setprecision(6) << i << " " << lats_outer[i] << " " << lons_outer[i] << " " << values_outer[i] << " " << flags_outer[i] << std::endl;
                }
                std::cout << "indices_outer_inner - l i lats lons obs:" << std::endl;
                for(int l=0; l<p_inner; l++) {
                    int i = indices_outer_inner[l];
                    std::cout << std::setprecision(6) << l << " " << i << " " << lats_outer[i] << " " << lons_outer[i] << " " << values_outer[i] << " " << flags_outer[i] << std::endl;
                }
                std::cout << "indices_outer_test - l m lats lons obs:" << std::endl;
                for(int m=0; m<p_test; m++) {
                    int i = indices_outer_test[m];
                    std::cout << std::setprecision(6) << m << " " << i << " " << lats_outer[i] << " " << lons_outer[i] << " " << values_outer[i] << " " << flags_outer[i] << std::endl;
                }
                std::cout << "indices_inner_test - l m lats lons obs:" << std::endl;
                for(int m=0; m<p_test; m++) {
                    int l = indices_inner_test[m];
                    int i = indices_outer_inner[l];
                    std::cout << std::setprecision(6) << m << " " << l << " " << lats_outer[i] << " " << lons_outer[i] << " " << values_outer[i] << " " << flags_outer[i] << std::endl;
                }
            }


            // -~- Compute the background 
            vec bvalues_outer = titanlibcustom::background(elevs_outer, values_outer, num_min_prof, min_elev_diff, value_minp, value_maxp, background_elab_type, background_values, indices_global_outer, debug);


            // *******************************************************************************************************************************
            // *******************************************************************************************************************************
            // *******************************************************************************************************************************
            // MN Edit: Print out background values
            if(debug){
                std::cout << "MN - Printing out the background values -- iteration: " << iteration << " -- observation: " << curr << std::endl;
                for(int i=0; i<bvalues_outer.size(); i++) {
                    std::cout << bvalues_outer[i] << std::endl;
                }
            }
            // *******************************************************************************************************************************
            // *******************************************************************************************************************************
            // *******************************************************************************************************************************

            if(debug) std::cout << "... background ok ..." << std::endl;

            // -~- If deviations between backgrounds and observations are small then flag = 0

            /* if observations and background are almost identical (within the range of valid estimates),
               then take a shortcut and flag = 0 for all observations in the inner circle */

            bool small_innov = true;
            for(int m=0; m<p_test; m++) {
                int j = indices_outer_test[m];
                if(bvalues_outer[j] < minv_outer[j] || bvalues_outer[j] > maxv_outer[j]) 
                    small_innov = false;
            }

            if (small_innov) {
                for(int m=0; m<p_test; m++) {
                    int g = indices_global_test[m];
                    flags[g] = 0;
                    if (debug) std::cout << " small_innov - index " << g << std::endl;
                }
                continue;
            }
 
            // -~- FGT on a selection of observations to test in the inner circle
            bool set_flag0 = iteration > 0;
            res = fgt_core( lats_outer, lons_outer, values_outer, bvalues_outer, sigmab_outer, value_minp, value_maxp, mina_outer, maxa_outer, minv_outer, maxv_outer, tpos_outer, tneg_outer, indices_global_test, indices_outer_inner, indices_outer_test, indices_inner_test, debug, na, set_flag0, thrown_out, scores, flags);
            // thrown_out is reset every new iteration

            // problems during the matrix inversion
            if ( !res) {
                flags[curr] = 100;
                if(debug) std::cout << " oi - flags=100 - index " << curr << std::endl; 
                continue;
            }

            count_tests++; // one more tests

        }  // end loop over observations

        std::cout << "FGT loop - Removing " << thrown_out << " observations. Number of tests " << count_tests << std::endl;
        double e_time0 = titanlibcustom::clock();
        std::cout << e_time0 - s_time0 << " secs" << std::endl;
        if(thrown_out == 0) {
            if ( iteration == 0) set_all_good = true;
            if(iteration + 1 < num_iterations)
                std::cout << "Stopping early after " << iteration + 1<< " iterations" << std::endl;
            break;
        }
    
    } // end of FGT iterations

    /* enter if the first iteration has been completed without flagging any bad observation
       The loop flags all the observations as good ones */
    if ( set_all_good) {
        int count_zero = 0;
        for(int curr=0; curr < p; curr++) {
          if ( flags[curr] == na && obs_test[curr] == 1) {
            flags[curr] = 0;
            count_zero++;
          }
        }
        if(debug) std::cout << "Remaining " << count_zero << " observations set to \"good\" " << std::endl;
    }

    /* check on observations missing QC flags */

    if(debug) std::cout << " +++++ check on observations missing QC flags ++++++++++++++++" << std::endl;
    double s_time0 = titanlibcustom::clock();
        
    // reset this number each loop (this is for breaking if we don't throw anything new out)
    int thrown_out = 0; 
        
    // diagnostic. count the number of times the test is performed
    int count_tests = 0;

    // loop over observations
    for(int curr=0; curr < p; curr++) {

        if(debug) std::cout << "===> curr " << curr << " ===============" << std::endl;

        // -~- Can the observation be used as a centroid for the definition of inner and outer circles?

        if(obs_test[curr] != 1 || flags[curr] >= 0) {
            if(debug) std::cout << "..skip " << curr << std::endl;
            continue;
        }

        // -~- Define outer/inner circles and which are the observations to test

        // get all neighbours that are close enough (inside outer circle)
        vec distances;
        ivec indices_global_outer_guess = tree.get_neighbours_with_distance(lats[curr], lons[curr], outer_radius, distances, true);
        if(indices_global_outer_guess.size() > num_max_outer) {

            int N = indices_global_outer_guess.size();
            std::vector<std::pair<float,int> > pairs(N);
            for(int i = 0; i < indices_global_outer_guess.size(); i++) {
                pairs[i] = std::pair<float, int>(distances[i], indices_global_outer_guess[i]);
            }
            std::sort(pairs.begin(), pairs.end(), titanlibcustom::sort_pair_first<float,int>());
            distances.clear();
            indices_global_outer_guess.clear();
            distances.resize(num_max_outer);
            indices_global_outer_guess.resize(num_max_outer);
            for(int i = 0; i < num_max_outer; i++) {
                distances[i] = pairs[i].first;
                indices_global_outer_guess[i] = pairs[i].second;
            }
        }

        // set all the indices linking the different levels
        // note that the only observation to test is curr
        ivec indices_global_outer, indices_global_test, indices_outer_inner, indices_outer_test, indices_inner_test;
        bool res = titanlibcustom::set_indices( indices_global_outer_guess, obs_test, flags, distances, inner_radius, curr, indices_global_outer, indices_global_test, indices_outer_inner, indices_outer_test, indices_inner_test);
            
        int p_outer = indices_global_outer.size();
        int p_inner = indices_outer_inner.size();
        int p_test = indices_global_test.size();
        if(debug) std::cout << "p_outer inner test " << p_outer << " " << p_inner << " " << p_test << std::endl;
            
        // -~- Check if there are enough observations in the outer/inner circles

        if(p_outer < num_min_outer) {
            flags[curr] = 12;
            if(debug) std::cout << "@@isolated (outer) " << curr << std::endl;
            continue; 
        }

        if( p_inner < 2) {
            flags[curr] = 11;
            if(debug) std::cout << "@@isolated (inner) " << curr << std::endl;
            continue;
        }
            
        // -~- Vectors on the outer circle

        vec lons_outer   = titanlibcustom::subset(lons, indices_global_outer);
        vec elevs_outer  = titanlibcustom::subset(elevs, indices_global_outer);
        vec lats_outer   = titanlibcustom::subset(lats, indices_global_outer);
        vec sigmab_outer = titanlibcustom::subset(sigma_b, indices_global_outer);
        vec values_outer = titanlibcustom::subset(values, indices_global_outer);
        vec tpos_outer   = titanlibcustom::subset(tpos, indices_global_outer);
        vec tneg_outer   = titanlibcustom::subset(tneg, indices_global_outer);
        vec mina_outer   = titanlibcustom::subset(values_mina, indices_global_outer);
        vec minv_outer   = titanlibcustom::subset(values_minv, indices_global_outer);
        vec maxa_outer   = titanlibcustom::subset(values_maxa, indices_global_outer);
        vec maxv_outer   = titanlibcustom::subset(values_maxv, indices_global_outer);


        // -~- Compute the background 
        vec bvalues_outer = titanlibcustom::background(elevs_outer, values_outer, num_min_prof, min_elev_diff, value_minp, value_maxp, background_elab_type, background_values, indices_global_outer, debug);

        if(debug) std::cout << "... background ok ..." << std::endl;

        // -~- If deviations between backgrounds and observations are small then flag = 0

        /* if observations and background are almost identical (within the range of valid estimates),
           then take a shortcut and flag = 0 for all observations in the inner circle */

        int j = indices_outer_test[0];
        if(bvalues_outer[j] > minv_outer[j] && bvalues_outer[j] < maxv_outer[j]) { 
            int g = indices_global_test[0];
            flags[g] = 0;
            if (debug) std::cout << " small_innov - index " << g << std::endl;
            continue;
        }

        // -~- FGT on a selection of observations to test in the inner circle

        res = fgt_core( lats_outer, lons_outer, values_outer, bvalues_outer, sigmab_outer, value_minp, value_maxp, mina_outer, maxa_outer, minv_outer, maxv_outer, tpos_outer, tneg_outer, indices_global_test, indices_outer_inner, indices_outer_test, indices_inner_test, debug, na, true, thrown_out, scores, flags);

        // problems during the matrix inversion
        if ( !res) {
            flags[curr] = 100;
            if(debug) std::cout << " oi - flags=100 - index " << curr << std::endl; 
            continue;
        }

        count_tests++; // one more tests

    }  // end loop over observations

    std::cout << "QC missing - Removing " << thrown_out << " observations. Number of tests " << count_tests << std::endl;
    double e_time0 = titanlibcustom::clock();
    std::cout << e_time0 - s_time0 << " secs" << std::endl;

    /* final check on the bad observations
       it may happen that a good observation is flagged as a bad one because of the order
       the FGT has been done and the uncertainty made in the estimation of the background 
       (remember that bad observations may have been used to get the background).
       This final check is made only on bad observations and uses only good observations. */

    if(debug) std::cout << " +++++ final check on the bad observations ++++++++++++++++" << std::endl;
    s_time0 = titanlibcustom::clock();
        
    // reset this number each loop (this is for breaking if we don't throw anything new out)
    thrown_out = 0; 
        
    // diagnostic. count the number of times the test is performed
    count_tests = 0;

    // loop over observations
    for(int curr=0; curr < p; curr++) {

        if(debug) std::cout << "===> curr " << curr << " ===============" << std::endl;

        // -~- Can the observation be used as a centroid for the definition of inner and outer circles?
        if(obs_test[curr] != 1 || flags[curr] != 1 || values[curr] < value_minp || values[curr] > value_maxp) {
            if(debug) std::cout << "..skip " << curr << std::endl;
            continue;
        }

        // -~- Define outer/inner circles and which are the observations to test

        // get all neighbours that are close enough (inside outer circle)
        vec distances;
        ivec indices_global_outer_guess = tree.get_neighbours_with_distance(lats[curr], lons[curr], outer_radius, distances, true);
        if(indices_global_outer_guess.size() > num_max_outer) {

            int N = indices_global_outer_guess.size();
            std::vector<std::pair<float,int> > pairs(N);
            for(int i = 0; i < indices_global_outer_guess.size(); i++) {
                pairs[i] = std::pair<float, int>(distances[i], indices_global_outer_guess[i]);
            }
            std::sort(pairs.begin(), pairs.end(), titanlibcustom::sort_pair_first<float,int>());
            distances.clear();
            indices_global_outer_guess.clear();
            distances.resize(num_max_outer);
            indices_global_outer_guess.resize(num_max_outer);
            for(int i = 0; i < num_max_outer; i++) {
                distances[i] = pairs[i].first;
                indices_global_outer_guess[i] = pairs[i].second;
            }
        }

        // set all the indices linking the different levels

        ivec indices_global_outer, indices_global_test, indices_outer_inner, indices_outer_test, indices_inner_test;
        bool res = titanlibcustom::set_indices(indices_global_outer_guess, obs_test, flags, distances, inner_radius, curr, indices_global_outer, indices_global_test, indices_outer_inner, indices_outer_test, indices_inner_test);
            
        int p_outer = indices_global_outer.size();
        int p_inner = indices_outer_inner.size();
        int p_test = indices_global_test.size();
        if(debug) std::cout << "p_outer inner test " << p_outer << " " << p_inner << " " << p_test << std::endl;
            
        // -~- Decide if there are enough observations in the outer/inner circles

        if(p_outer < num_min_outer) {
            flags[curr] = 12;
            if(debug) std::cout << "@@isolated (outer) " << curr << std::endl;
            continue; 
        }

        if( p_inner < 2) {
            flags[curr] = 11;
            if(debug) std::cout << "@@isolated (inner) " << curr << std::endl;
            continue;
        }
            
        // -~- Set vectors on the outer circle

        vec lons_outer   = titanlibcustom::subset(lons, indices_global_outer);
        vec elevs_outer  = titanlibcustom::subset(elevs, indices_global_outer);
        vec lats_outer   = titanlibcustom::subset(lats, indices_global_outer);
        vec sigmab_outer = titanlibcustom::subset(sigma_b, indices_global_outer);
        vec values_outer = titanlibcustom::subset(values, indices_global_outer);
        vec tpos_outer   = titanlibcustom::subset(tpos, indices_global_outer);
        vec tneg_outer   = titanlibcustom::subset(tneg, indices_global_outer);
        vec mina_outer   = titanlibcustom::subset(values_mina, indices_global_outer);
        vec minv_outer   = titanlibcustom::subset(values_minv, indices_global_outer);
        vec maxa_outer   = titanlibcustom::subset(values_maxa, indices_global_outer);
        vec maxv_outer   = titanlibcustom::subset(values_maxv, indices_global_outer);

        // Debug for indices
        if(debug) {
            ivec flags_outer;
            flags_outer.resize( p_outer, na);
            for(int i=0; i < p_outer; i++) {
               int g = indices_global_outer[i];
               flags_outer[i] = flags[g];
            }
            std::cout << "curr lats lons obs:" << std::endl;
            std::cout << curr << " " << lats[curr] << " " << lons[curr] << " " << values[curr] << std::endl;
            std::cout << "indices_global_outer - i g lats lons obs:" << std::endl;
            for(int i=0; i<p_outer; i++) {
                int g = indices_global_outer[i];
                std::cout << std::setprecision(6) << i << " " << g << " " << lats[g] << " " << lons[g] << " " << values[g] << " " << flags[g] << std::endl;
            }
            std::cout << "outer - i lats lons obs:" << std::endl;
            for(int i=0; i<p_outer; i++) {
                std::cout << std::setprecision(6) << i << " " << lats_outer[i] << " " << lons_outer[i] << " " << values_outer[i] << " " << flags_outer[i] << std::endl;
            }
            std::cout << "indices_outer_inner - l i lats lons obs:" << std::endl;
            for(int l=0; l<p_inner; l++) {
                int i = indices_outer_inner[l];
                std::cout << std::setprecision(6) << l << " " << i << " " << lats_outer[i] << " " << lons_outer[i] << " " << values_outer[i] << " " << flags_outer[i] << std::endl;
            }
            std::cout << "indices_outer_test - l m lats lons obs:" << std::endl;
            for(int m=0; m<p_test; m++) {
                int i = indices_outer_test[m];
                std::cout << std::setprecision(6) << m << " " << i << " " << lats_outer[i] << " " << lons_outer[i] << " " << values_outer[i] << " " << flags_outer[i] << std::endl;
            }
            std::cout << "indices_inner_test - l m lats lons obs:" << std::endl;
            for(int m=0; m<p_test; m++) {
                int l = indices_inner_test[m];
                int i = indices_outer_inner[l];
                std::cout << std::setprecision(6) << m << " " << l << " " << lats_outer[i] << " " << lons_outer[i] << " " << values_outer[i] << " " << flags_outer[i] << std::endl;
            }
        }

        // -~- Compute the background 
        vec bvalues_outer = titanlibcustom::background(elevs_outer, values_outer, num_min_prof, min_elev_diff, value_minp, value_maxp, background_elab_type, background_values, indices_global_outer, debug);

        if(debug) { 
            int j = indices_outer_test[0];
            for(int i=0; i < p_outer; i++) {
                std::cout << std::setprecision(6) << "values backg minv maxv " << values_outer[i] << " " << bvalues_outer[i] << " " << minv_outer[i] << " " << maxv_outer[i] << " " << std::endl; 
            }
            std::cout << "j " << j << std::endl; 
        }
        if(debug) std::cout << "... background ok ..." << std::endl;

        // -~- If deviations between backgrounds and observations are small then flag = 0

        /* if observations and background are almost identical (within the range of valid estimates),
           then take a shortcut and flag = 0 for all observations in the inner circle */

        int j = indices_outer_test[0];
        if(bvalues_outer[j] > minv_outer[j] && bvalues_outer[j] < maxv_outer[j]) { 
            int g = indices_global_test[0];
            flags[g] = 0;
            if (debug) std::cout << " small_innov - index " << g << std::endl;
            continue;
        }

        // -~- FGT within the inner circle
        // note that the only observation to test is curr
        res = fgt_core( lats_outer, lons_outer, values_outer, bvalues_outer, sigmab_outer, value_minp, value_maxp, mina_outer, maxa_outer, minv_outer, maxv_outer, tpos_outer, tneg_outer, indices_global_test, indices_outer_inner, indices_outer_test, indices_inner_test, debug, na, true, thrown_out, scores, flags);

        // problems during the matrix inversion
        if ( !res) {
            flags[curr] = 100;
            if(debug) std::cout << " oi - flags=100 - index " << curr << std::endl; 
            continue;
        }

        count_tests++; // one more tests

    }  // end loop over observations

    std::cout << "Re-check bad obs - Removing " << thrown_out << " observations. Number of tests" << count_tests << std::endl;
    e_time0 = titanlibcustom::clock();
    std::cout << e_time0 - s_time0 << " secs" << std::endl;

    //
    if(debug) {
        int count_bad = 0;
        int count_good = 0;
        int count_missing = 0;
        int count_iso_inner = 0;
        int count_iso_outer = 0;
        int count_fail_matinv = 0;
        int count_impossible = 0;
        // loop over observations
        for(int curr=0; curr < p; curr++) {
            if( flags[curr] == 1) {
                count_bad++;
            } else if( flags[curr] == 0) {
                count_good++;
            } else if ( flags[curr] == na) {
                count_missing++;
            } else if( flags[curr] == 11) {
                count_iso_inner++;
            } else if( flags[curr] == 12) {
                count_iso_outer++;
            } else if( flags[curr] == 100) {
                count_fail_matinv++;
            } else {
                count_impossible++;
            }
        }
        std::cout << std::setprecision(3) << "summary - # TOT good bad missing isolated(inner) isolated(outer): " << p << " " << count_good << " " << count_bad << " " << count_missing << " " << count_iso_inner << " " << count_iso_outer << std::endl; 
        if(count_fail_matinv > 0) 
            std::cout << std::setprecision(3) << "!!!! failure in matrix inversion: " << count_fail_matinv << std::endl; 
        if(count_impossible > 0) 
            std::cout << std::setprecision(3) << "!!!! unknown flag: " << count_impossible << std::endl; 
    }
    //
    
    std::cout << ">> Total Time " << e_time0 - s_time << "secs" << std::endl;

    //
    return flags;
}
// end FGT //


//=============================================================================

//----------------------------------------------------------------------------//
// FGT WITHIN THE INNER CIRCLE

// Smallets unit for the first-guess test
bool fgt_core( const vec& lats, 
               const vec& lons, 
               const vec& yo, 
               const vec& yb, 
               const vec& sigma_b, 
               float minp,
               float maxp,
               const vec& mina, 
               const vec& maxa, 
               const vec& minv, 
               const vec& maxv, 
               const vec& tpos, 
               const vec& tneg, 
               const ivec& indices_global_test, 
               const ivec& indices_outer_inner, 
               const ivec& indices_outer_test, 
               const ivec& indices_inner_test, 
               bool debug, 
               float na, 
               bool set_flag0, 
               int& thrown_out,
               vec& scores,
               ivec& flags) {
/*-----------------------------------------------------------------------------
 First guess test for a selection of observations within the inner circle.

 The inner circle is defined with respect to a centroid observation.

 Each observation to test is compared against an estimated value (first guess)
 obtained considering the neighbouring observations.
 
 FGT scores:
 chi = abs( observation - background) / uncertainty of the background 
 The background uncertainty is optional, when not provided it is set to 1.
 z = (chi - mu) / (sigma + sigma_mu)
   mu = median( chi(*))
   sigma = inter-quartile range of chi(*)
   sigma_mu = sigma / sqrt(n(*))
 (*) statistics based on the n observations that are within the inner circle
 and with a background within the range of admissible values.

 FGT fails when:
 z > tpos, yo >= yb
 z > tneg, yo <  yb

 Our prior knowledge includes the definition of two confidence levels for the 
 background:
  - valid values, range [minv,maxv]
      background is so close that it validates the observation
  - admissible values, range [mina,maxa] 
      if background is outside this range, then the observation is bad 
  Note that the statistics defining the spread used in the FGT is based on 
  deviations between observations and estimates that lie within the range of 
  admissible values only.

 Decision tree:

 ALL observations are set to bad when:
  the backgrounds all lie outside the range of admissible values 

 ONE observation is deemed as bad when:
  fails the FGT with the largest z-score 
  AND 
  background is outside the range of valid values

 ALL observations are set to good when:
  the user allows the procedure to flag them as good ones (set_flag0 is true)
  AND
   the backgrounds all lie within the range of valid values
   OR
   none of the observations fail the FGT test

  Output values:

   flags, global vector. return updated quality control flag values (1 = bad, 0 = good)
   
   scores, global vector. return updated scores. Note that only for the case of 
                          a single bad observation we store the z-score
   
   thrown_out, integer. update the number of observations flagged as bad ones 

-----------------------------------------------------------------------------*/

    using namespace boost::numeric::ublas;

    // init
//    int p_outer = yo.size();
    int p_inner = indices_outer_inner.size();
    int p_test = indices_global_test.size();

    // obtain chi = sqrt( analysis_residual * cvanalysis_residual) 
    // two chi vectors are obtained:
    // 1. chi_inner, chi for all observations in the inner circle
    // 2. chi_stat, chi values used to extract summary stastics 
    //    observations in the inner circle AND yav has an admissible value   
    vec chi_stat;
    vec chi_stat_alt;
    chi_stat.reserve(p_inner);
    chi_stat_alt.reserve(p_inner);
    boost::numeric::ublas::vector<float> chi_inner(p_inner);
    boost::numeric::ublas::vector<float> chi_inner_alt(p_inner);
    if(debug) std::cout << "collect statistics - i l lon lat yo yb sigma_b chi chi_alt" << std::endl;
    for(int l=0; l<p_inner; l++) {
        int i = indices_outer_inner[l];
        chi_inner(l) = fabs( yo[i] - yb[i]) / sigma_b[i]; 
        chi_inner_alt(l) = maxv[i] - minv[i];
        if ( yb[i] >= mina[i] && yb[i] <= maxa[i]) { 
            chi_stat.push_back(chi_inner(l));
            chi_stat_alt.push_back(chi_inner_alt(l));
        }
        if(debug) std::cout << std::setprecision(6) << i << " " << l << " " << lons[i] << " " << lats[i] << " " << yo[i] << " " << yb[i] << " " << sigma_b[i] << " " << chi_inner(l) << " " << chi_inner_alt(l) << std::endl;
    }
    
    // chi_stat is empty (yav all outside the range of admissible values), set all flags to 1
    if(chi_stat.size() == 0) {
        for(int m=0; m<p_test; m++) {
            int g = indices_global_test[m];
            flags[g] = 1;
            thrown_out++;
            if(debug) std::cout << std::setprecision(3) << "chi_stat empty - flag as bad " << g << std::endl;
        }
        return true;       
    }

    // chi summary statistics
    float mu = titanlibcustom::compute_quantile( 0.5, chi_stat);
    float sigma = titanlibcustom::compute_quantile( 0.75, chi_stat) - titanlibcustom::compute_quantile( 0.25, chi_stat);
    float sigma_alt = titanlibcustom::compute_quantile( 0.75, chi_stat_alt) - titanlibcustom::compute_quantile( 0.25, chi_stat_alt);


    if (sigma_alt > sigma) sigma = sigma_alt;
 
    if(sigma == 0) {
        if(debug) std::cout << "sigma = 0, not possible to perform the spatial check " << std::endl;
        return true;       
    }
    float sigma_mu = sigma / std::sqrt( chi_stat.size());

    if(debug) {
        std::cout << " chi_stat, dim = " << chi_stat.size() << std::endl;
        for(int l=0; l<chi_stat.size(); l++) 
            std::cout << std::setprecision(6) << chi_stat[l] << std::endl;
        std::cout << std::setprecision(6) << " mu sigma sigma_mu " << mu << " " << sigma << " " << sigma_mu << std::endl;
    }
   
    // z = ( chi - mu ) / ( sigma + sigma_mu )
    float zmx = -10000;
    int mmx = -1;
    if(debug) std::cout << "z_test - i l lons lats yo yb chi z " << std::endl;
    for(int m=0; m<p_test; m++) {
        int i = indices_outer_test[m];
        int l = indices_inner_test[m];
        float z = (chi_inner(l) - mu) / (sigma + sigma_mu);
        if (z > zmx && (yb[i] < minv[i] || yb[i] > maxv[i])) {
          zmx = z;
          mmx = m;
        } 
        if(debug) std::cout << std::setprecision(6) << i << " " << l << " " << lons[i] << " " << lats[i] << " " << yo[i] << " " << yb[i] << " " << chi_inner(l) << " " << z << std::endl;
    }

    // Decision making
    float thr = 0;
    int i = 0;
    int l = 0;
    if ( mmx >= 0) { 
        i = indices_outer_test[mmx];
        l = indices_inner_test[mmx];
        // set the threshold
        thr = tneg[i];
        if ( (yo[i] - yb[i]) >= 0) thr = tpos[i];
    }
    // if the largest z is larger than the threshold, then flag IT as bad ...
    if ( zmx > thr) {
        int g = indices_global_test[mmx];
        scores[g] = zmx;
        flags[g] = 1;
        thrown_out++;
        if(debug) std::cout << std::setprecision(6) << "FGT failed - flag as bad - index yo yb chi z " << g << " " << yo[i] << " " << yb[i] << " " << chi_inner(l) << " " << zmx << std::endl;
    // ... BUT if the largest z is smaller than the threshold (and we are confindet in flagging), then flag ALL as good
    } else if (set_flag0) {
        for(int m=0; m<p_test; m++) {
            int g = indices_global_test[m];
            flags[g] = 0;
            if(debug) std::cout << "all set to good - g " << g << std::endl;
        }
    }
    
    // normal exit
    return true;       
}
