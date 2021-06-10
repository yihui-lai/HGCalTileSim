/**
 * @file MCHelper.cc
 * @author Yi-Mu "Enoch" Chen
 * @brief Implementing the functions required for MC topology crawling.
 */
#include "DataFormats/Candidate/interface/Candidate.h"

#include <queue>

namespace usr {

/**
 * @brief returning pointer to the particle in the parent topology with a
 *        specified flavour.
 *
 * For an input particle x, it will traverse the topology towards the parent
 * side, until
 * 1. the parent particle has the target flavour F, in which case it returns
 *    the pointer to the parent particle. Note that the absolute value will
 *    be taken for the for the flavour.
 * 2. The parent particle has a flavour different to x, in which case a null
 *    pointer is return.
 *
 *  This function help find direct parent of a particle while avoid the
 *  possibility of final state radiations affecting the topology parsing. For
 *  example, the following topology is actually rather common in CMSSW truth
 *  topologies:
 *
 * \f[
 *      F->(x+g)->(x'+gg)->(x''+ggg)
 * \f]
 *  The three function calls GetDirectMother( x, F ),
 *  GetDirectMother( x', F ), and GetDirectMother( x'', F ) should all return
 *  the same pointer F
 *
 *  Topologies in MiniAOD also have some times have the strange topologies like:
 *  x -> x' -> x'' (with no additional radiation particles) which could be
 *  avoided using this function.
 */
const reco::Candidate*
GetDirectMother( const reco::Candidate* x, int target_ID )
{
  std::queue<const reco::Candidate*> bfs_queue;
  bfs_queue.push( x );

  while( !bfs_queue.empty() ){
    const reco::Candidate* temp = bfs_queue.front();
    bfs_queue.pop();
    if( abs( temp->pdgId() ) == abs( target_ID ) ){
      return temp;
    }

    for( unsigned i = 0; i < temp->numberOfMothers(); ++i ){
      if( temp->mother( i )->pdgId() == temp->pdgId() ){
        bfs_queue.push( temp->mother( i ) );
      } else if( abs( temp->mother( i )->pdgId() ) == abs( target_ID ) ){
        return temp->mother( i );
      }
    }
  }

  return nullptr;
}

/**
 * @brief Essentially the same as the GetMother functions, except this time the
 *        traversal directions is towars the children side.
 *
 * It will follow the same flavour decaying chain until there is a flavour
 * change. A major caveat of this functions would be if the decay topology
 * you are interested in is on the otherside of the hard radiation photon,
 * gluon or Z boson. In these cases, the user should be careful that this
 * function is returning what is desired.
 */
const reco::Candidate*
GetDaughter( const reco::Candidate* x, int target_ID )
{
  std::queue<const reco::Candidate*> bfs_queue;
  bfs_queue.push( x );

  while( !bfs_queue.empty() ){
    const reco::Candidate* temp = bfs_queue.front();
    bfs_queue.pop();
    if( abs( temp->pdgId() ) == abs( target_ID ) ){
      // Moving to bottom of single decay chain
      while( temp->numberOfDaughters() == 1 ){
        temp = temp->daughter( 0 );
      }

      return temp;
    }

    for( unsigned i = 0; i < temp->numberOfDaughters(); ++i ){
      bfs_queue.push( temp->daughter( i ) );
    }
  }

  return nullptr;
}

}/* usr */
