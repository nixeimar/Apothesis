#include "abstract_process.h"
#include "factory_process.h"

AbstractProcess::AbstractProcess( const string& name ){
    FactoryProcess::registerThis( name, this);
}
