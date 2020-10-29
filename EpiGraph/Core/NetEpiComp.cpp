//
// Created by roryh on 28/10/2020.
//

#include <EpiGraph/Core/NetEpiComp.hpp>

namespace EpiGraph {

    template<>
    auto sixrd_meta_pop_ode(const NetEpiComp<1> &model) -> typename NetEpiComp<1>::state_type {
        return net_SIXRD_ode(model.state(), model.coupling(), model.params().transpose());
    }

    template<>
    auto sixrd_meta_pop_ode(const NetEpiComp<Eigen::Dynamic> &model) -> typename NetEpiComp<Eigen::Dynamic>::state_type {
        return net_SIXRD_ode_inhom(model.state(), model.coupling(), model.params());
    }

    // both are instantiated above
    //template<> class NetEpiComp<1>;
    //template<> class NetEpiComp<Eigen::Dynamic>;

}