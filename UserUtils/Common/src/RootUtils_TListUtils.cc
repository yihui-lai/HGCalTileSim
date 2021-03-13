#include "TList.h"

#include <stack>

namespace usr {

/**
 * @brief Setting the alls the `target` instances in a TList object to be placed
 * before the first instance of `before` in the TList.
 *
 * This function is required since the public member functions of TList doesn't
 * allow easy manipulation of the accompanied options string. If `before` doesn't
 * exist in the TList, the function will do nothing and return false. Otherwise
 * the function will return true, regardless of whether anything about the TList
 * is actually changed by the function or not.
 *
 * The implementation of the function is basically dumping all contents after the
 * first `before` instance into a stack, appended the required `target`
 * instances, and reinserting all of the other instances.
 */
bool
MoveObjectToBefore( TList&         list,
                    const TObject& target,
                    const TObject& before )
{
  struct objopt
  {
    TObject*    obj;
    std::string opt;
  };

  // Finding the link to where the targets should be moved.
  TObjLink* beforelink = list.FirstLink();

  while( beforelink && beforelink->GetObject() != &before ){
    beforelink = beforelink->Next();
  }

  if( !beforelink ){ return false; }
  beforelink = beforelink->Prev();

  // Flushing all the stuff behind into a target position.
  std::stack<objopt> otherlinks;
  std::stack<objopt> targetlinks;

  TObjLink* link = list.LastLink();

  while( link != beforelink ){
    const objopt temp = { link->GetObject(), link->GetOption() };
    if( link->GetObject() == &target ){
      targetlinks.push( temp );
    } else {
      otherlinks.push( temp );
    }
    link = link->Prev();
    list.RemoveLast();
  }

  // Inserting targets
  while( !targetlinks.empty() ){
    list.Add(
      targetlinks.top().obj,
      targetlinks.top().opt.c_str() );
    targetlinks.pop();
  }

  // Inserting the rest of the stuff.
  while( !otherlinks.empty() ){
    list.Add(
      otherlinks.top().obj,
      otherlinks.top().opt.c_str() );
    otherlinks.pop();
  }

  return true;

}

} /* namespace usr */