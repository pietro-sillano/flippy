#ifndef FLIPPY_MONTECARLOUPDATER_HPP
#define FLIPPY_MONTECARLOUPDATER_HPP
#include "custom_concepts.hpp"
#include <random>
#include "Nodes.hpp"
#include "Triangulation.hpp"

namespace fp {

/**
 * @brief A helper class for updating of the triangulation, using [Metropolis–Hastings algorithm](https://en.wikipedia.org/wiki/Metropolis%E2%80%93Hastings_algorithm).
 *
 * @tparam Real @RealStub
 * @tparam Index @IndexStub
 * @tparam EnergyFunctionParameters
 * @tparam RandomNumberEngine
 * @tparam triangulation_type
 */
template<floating_point_number Real, indexing_number Index, typename EnergyFunctionParameters, typename RandomNumberEngine, TriangulationType triangulation_type>
class MonteCarloUpdater
{
private:
    static constexpr Real max_float = 3.40282347e+38;
    Real e_old{}, e_new{}, e_diff{};
    fp::Triangulation<Real, Index, triangulation_type>& triangulation;
    EnergyFunctionParameters const& prms;
    std::function<Real(fp::Node<Real, Index> const&, fp::Triangulation<Real, Index, triangulation_type> const&, EnergyFunctionParameters const&)> energy_function;
    RandomNumberEngine& rng;
    std::uniform_real_distribution<Real> unif_distr_on_01;
    Real kBT_{1};
    Real min_bond_length_square{0.}, max_bond_length_square{max_float};
    unsigned long move_attempt{0}, bond_length_move_rejection{0},move_back{0};
    unsigned long flip_attempt{0}, bond_length_flip_rejection{0}, flip_back{0};

public:
    MonteCarloUpdater(fp::Triangulation<Real, Index, triangulation_type>& triangulation_inp,
                      EnergyFunctionParameters const& prms_inp,
                      std::function<Real(fp::Node<Real, Index> const&, fp::Triangulation<Real, Index, triangulation_type> const&, EnergyFunctionParameters const&)> energy_function_inp,
                      RandomNumberEngine& rng_inp, Real min_bond_length, Real max_bond_length)
    :triangulation(triangulation_inp), prms(prms_inp), energy_function(energy_function_inp), rng(rng_inp),
    unif_distr_on_01(std::uniform_real_distribution<Real>(0, 1)),
    min_bond_length_square(min_bond_length*min_bond_length), max_bond_length_square(max_bond_length*max_bond_length)
    {}

    bool move_needs_undoing()
    {
        e_diff = e_old - e_new;
        if(kBT_>0){ //temperature can safely be put to 0, this will make the algorithm greedy
            return (e_diff<0) && (unif_distr_on_01(rng)>std::exp(e_diff/kBT_));
        }else{
            return (e_diff<0);
        }
    }

    bool new_neighbour_distances_are_between_min_and_max_length(fp::Node<Real, Index> const& node,
                                                                     fp::vec3<Real> const& displacement)

    {
        return (new_next_neighbour_distances_are_between_min_and_max_length(node, displacement)&&
            new_verlet_neighbour_distances_are_between_min_and_max_length(node, displacement));

    }

    bool new_next_neighbour_distances_are_between_min_and_max_length(fp::Node<Real, Index> const& node,
                                                                     fp::vec3<Real> const& displacement)

    {
        /**
         * iterate through the next neighbour distances of a node which is a collection of
         * 3d vectors, i.e. nn_dist is of type fp::vec3. Then check if the displacement would make any of the
         * nn_distance vectors longer than the allowed max length.
         */

        Real distance_square_new, distance_square_old;
        for (auto const& nn_dist: node.nn_distances){
            distance_square_new=(nn_dist - displacement).norm_square();
            distance_square_old=nn_dist.norm_square();
            if ((distance_square_new>max_bond_length_square) && (distance_square_old < max_bond_length_square)) {
                return false;
            }
            if ((distance_square_old>min_bond_length_square) && (distance_square_new<min_bond_length_square)) {
                return false;
            }
        }

        return true;
    }

    bool new_verlet_neighbour_distances_are_between_min_and_max_length(fp::Node<Real, Index> const& node,
                                                                     fp::vec3<Real> const& displacement)

    {
        /**
         * iterate through the next neighbour distances of a node which is a collection of
         * 3d vectors, i.e. nn_dist is of type fp::vec3. Then check if the displacement would make any of the
         * nn_distance vectors longer than the allowed max length.
         */
        Real distance_square_new, distance_square_old;
        for (auto const& verlet_neighbour_id: node.verlet_list)
        {
            distance_square_new=(triangulation[verlet_neighbour_id].pos - node.pos - displacement).norm_square();
            distance_square_old=(triangulation[verlet_neighbour_id].pos - node.pos).norm_square();
            if ((distance_square_new<min_bond_length_square)&&(distance_square_old>min_bond_length_square)) { return false; }
        }
        return true;
    }

    void move_MC_updater(fp::Node<Real, Index> const& node, fp::vec3<Real> const& displacement)
    {
        ++move_attempt;
        if (new_neighbour_distances_are_between_min_and_max_length(node, displacement)) {
            e_old = energy_function(node, triangulation, prms);
            triangulation.move_node(node.id, displacement);
            e_new = energy_function(node, triangulation, prms);
            if (move_needs_undoing()) {triangulation.move_node(node.id, -displacement); ++move_back;}
        }else{++bond_length_move_rejection;}
    }

    void flip_MC_updater(fp::Node<Real, Index> const& node)
    {
        ++flip_attempt;
        e_old = energy_function(node, triangulation, prms);
        Index number_nn_ids = node.nn_ids.size();
        Index nn_id = node.nn_ids[std::uniform_int_distribution<Index>(0, number_nn_ids-1)(rng)];
        auto bfd = triangulation.flip_bond(node.id, nn_id, min_bond_length_square, max_bond_length_square);
        if (bfd.flipped) {
            e_new = energy_function(node, triangulation, prms);
            if (move_needs_undoing()) { triangulation.unflip_bond(node.id, nn_id, bfd); ++flip_back;}
        }else{++bond_length_flip_rejection;}
    }

    void flip_MC_updater(fp::Node<Real, Index> const& node, Index index_in_nn_ids)
    {
        /**
         *
         */
        ++flip_attempt;
        e_old = energy_function(node, triangulation, prms);
        Index number_nn_ids = node.nn_ids.size();
//        Index nn_id = index_in_nn_ids;//node.nn_ids[std::uniform_int_distribution<Index>(0, number_nn_ids-1)(rng)];
        auto bfd = triangulation.flip_bond(node.id, index_in_nn_ids, min_bond_length_square, max_bond_length_square);
        if (bfd.flipped) {
            e_new = energy_function(node, triangulation, prms);
            if (move_needs_undoing()) { triangulation.unflip_bond(node.id, index_in_nn_ids, bfd); ++flip_back;}
        }else{++bond_length_flip_rejection;}
    }


    void reset_kBT(Real kBT){
        /**
         * @param kBT input value of new temperature.
         * The internal private data member `kBT_` that contains the current value of temperature will be overwritten
         * with this input.
         */
        kBT_=kBT;
    }

    [[nodiscard]] Real kBT(){
    /**
     * Monte Carlo Updater requires a kBT value in teh calculation of the rejection probability of a random move.
     * @return the state of kBT value of the updater.
     * @see move_needs_undoing()
     */
        return kBT_;
    }
    [[nodiscard]] unsigned long move_attempt_count() const {
    /**
     * Every time a move is attempted, a private internal state variable `move_attempt` is incremented by move_MC_updater().
     * This variable can be used for diagnostics or statistical tracking but its state does not impact the function of the updater.
     * @return current state of `move_attempt`.
     * @see move_needs_undoing()
     */
        return move_attempt;
    }
    [[nodiscard]] unsigned long bond_length_move_rejection_count() const {
    /**
     * Moves that cause nodes to overlap with their Verlet list neighbours or move them too far away from any of their next neighbours, are rejected. Every time such rejection happens, a private internal state variable `bond_length_move_rejection` is incremented by move_MC_updater().
     * The specifics of this rejection criteria are calculated in the function new_neighbour_distances_are_between_min_and_max_length(fp::Node<Real, Index> const& node, fp::vec3<Real> const& displacement)
     * Every time a node move would lead that node to have a bond with one of its next neighbours which is longer than a specified maximal length (max_bond_length()), the move will be rejected.
     * This variable can be used for diagnostics or statistical tracking but its state does not impact the function of the updater.
     * @return current state of `bond_length_move_rejection`.
     * @see move_needs_undoing()
     */
        return bond_length_move_rejection;
    }
    [[nodiscard]] unsigned long move_back_count() const {
    /**
     * Every time a move is rejected, because the energy requirement was not satisfied, a private internal state variable `move_back` is incremented by move_MC_updater().
     * This variable does not track the rejections resulting from bond length restriction violation.
     * This variable can be used for diagnostics or statistical tracking but its state does not impact the function of the updater.
     * @return current state of `move_back`.
     * @see move_needs_undoing() bond_length_move_rejection_count()
     */
        return move_back;
    }
    [[nodiscard]] unsigned long flip_attempt_count() const {
    /**
     * Every time a flip is attempted, a private internal state variable `flip_attempt` is incremented by flip_MC_updater() and flip_MC_updater(fp::Node<Real, Index> const& node, Index index_in_nn_ids).
     * This variable can be used for diagnostics or statistical tracking but its state does not impact the function of the updater.
     * @return current state of `flip_attempt`.
     */
        return flip_attempt;
    }
    [[nodiscard]] unsigned long bond_length_flip_rejection_count() const {
    /**
     * If a flip would turn a valid bond into a bond that is too long, the flip is rejected, a private internal state variable `bond_length_flip_rejection` is incremented by flip_MC_updater() and flip_MC_updater(fp::Node<Real, Index> const& node, Index index_in_nn_ids).
     * The rejection of flips is handled by the Triangulation class itself and reported through the BondFlipData to the flip functions of the updater.
     * This variable can be used for diagnostics or statistical tracking but its state does not impact the function of the updater.
     * @return current state of `bond_length_flip_rejection`.
     * @see Triangulation::flip_bond(Index, Index, Real, Real) Triangulation::unflip_bond(Index, Index, BondFlipData<Index> const&)
     */
        return bond_length_flip_rejection;
    }
    [[nodiscard]] unsigned long flip_back_count() const {
    /**
     * Every time a flip is rejected, because the energy requirement was not satisfied, a private internal state variable `flip_back` is incremented by flip_MC_updater() and flip_MC_updater(fp::Node<Real, Index> const& node, Index index_in_nn_ids).
     * The rejection of flips is handled by the Triangulation class itself and reported through the BondFlipData to the flip functions of the updater.
     * This variable does not track the rejections resulting from bond length restriction violation.
     * This variable can be used for diagnostics or statistical tracking but its state does not impact the function of the updater.
     * @return current state of `flip_back`.
     * @see Triangulation:: unflip_bond(Index node_id, Index nn_id, BondFlipData<Index> const& common_nns) bond_length_flip_rejection_count()
     */
        return flip_back;
    }

};
}
#endif //FLIPPY_MONTECARLOUPDATER_HPP
