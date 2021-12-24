

#pragma once

#include <kitty/constructors.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/operations.hpp>
#include "../utils/node_map.hpp"
#include "miter.hpp"
#include "simulation.hpp"
#include <cmath>

namespace mockturtle
{

/* Statistics to be reported */
struct simulation_cec_stats{
  //Here we write  Split Variable, that is related to the simulation size



  uint32_t split_var{ 0 };


  //Number of simulation rounds

  uint32_t rounds{ 0 };
};

namespace detail{

   //start here---------------------------------------------------------------------

class split_var_renato{
public:

    unsigned num_variables;

    unsigned split_variable;

    uint64_t others_round;
    //Destructor of the class
  split_var_renato() = delete;

  //This is the contructor of the class
  split_var_renato( unsigned num_variables, unsigned split_variable, uint64_t others_round ) : num_variables( num_variables ), split_variable{ split_variable }, others_round{ others_round } {}

  kitty::dynamic_truth_table compute_constant( bool value ) const
  {
    kitty::dynamic_truth_table tt( split_variable );
    return value ? ~tt : tt;
  }

  kitty::dynamic_truth_table compute_pi( uint32_t index ) const
  {
    kitty::dynamic_truth_table tt( split_variable );
    if ( index < split_variable ) {                        
      kitty::create_nth_var( tt, index );
    }
    else   {
     
      bool value = ( others_round >> ( index - split_variable ) ) & 1;
      if ( !value ){
        tt = ~tt;


      }
    }

    return tt;
  }

  kitty::dynamic_truth_table compute_not( kitty::dynamic_truth_table const& value ) const
  {
    return ~value;
  }

private:
 
  
  //uint64_t others_round;
};
//finish here --------------------------------------------------
template<class Ntk>
class simulation_cec_impl
{
public:
  using pattern_t = unordered_node_map<kitty::dynamic_truth_table, Ntk>;

  using node = typename Ntk::node;
  using signal = typename Ntk::signal;

public:

  explicit simulation_cec_impl( Ntk& ntk, simulation_cec_stats& st )
      : _ntk( ntk ),
        _st( st )
  {
  }

  bool run()
  {
 
    // START HERE ----------------------
    unsigned n = _ntk.num_pis();
    unsigned split_variable;
    unsigned V = _ntk.size();

   
    if ( n <= 6 )
    {
      split_variable = n;
      //or also _st.split_variable
   }
    else {
     

      int m = 7;
      while ( m < n && ( 32 + ( 1 << ( ( m + 1 ) - 3 ) ) ) * V <= 1 << 29 )
      {
        m++;
      }
      split_variable = m;
      //or _st.split_variable
    }

    int rounds = 1 << ( n - split_variable );

    // I store them in the statistics struct
    _st.split_var = split_variable;
    _st.rounds = rounds;

    // This is the actual simulation...

    for ( uint64_t others_round = 0; others_round < rounds; others_round++ )
    { // We iterate over all the possible assignations of the remaining variables
      split_var_renato simul( _ntk.num_pis(), split_variable, others_round );
      const auto tts = simulate<kitty::dynamic_truth_table>( _ntk, simul );

      // const std::vector<kitty::dynamic_truth_table> tts
      for ( auto& po : tts )
      {
        if ( !kitty::is_const0( po ) )
        {
          return false;
        }
      }
    }

    return true;
  }
   

  //FINISH HERE ---------------------------


private:
  Ntk& _ntk;
  simulation_cec_stats& _st;
  
};

} 


template<class Ntk>
std::optional<bool> simulation_cec( Ntk const& ntk1, Ntk const& ntk2, simulation_cec_stats* pst = nullptr )
{
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( has_num_pis_v<Ntk>, "Ntk does not implement the num_pis method" );
  static_assert( has_size_v<Ntk>, "Ntk does not implement the size method" );
  static_assert( has_get_node_v<Ntk>, "Ntk does not implement the get_node method" );
  static_assert( has_foreach_pi_v<Ntk>, "Ntk does not implement the foreach_pi method" );
  static_assert( has_foreach_po_v<Ntk>, "Ntk does not implement the foreach_po method" );
  static_assert( has_foreach_node_v<Ntk>, "Ntk does not implement the foreach_node method" );
  static_assert( has_is_complemented_v<Ntk>, "Ntk does not implement the is_complemented method" );

  simulation_cec_stats st;

  bool result = false;

  if ( ntk1.num_pis() > 40 )
    return std::nullopt;

  auto ntk_miter = miter<Ntk>( ntk1, ntk2 );

  if ( ntk_miter.has_value() )
  {
    detail::simulation_cec_impl p( *ntk_miter, st );
    result = p.run();
  }

  if ( pst )
    *pst = st;

  return result;
}



} 
