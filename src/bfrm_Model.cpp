#include <RcppCommon.h>
#include <Model.h>
#include <Rcpp.h>

using namespace Rcpp ;

RCPP_MODULE(Model_module) {
	class_<Model>("Model")
	.constructor()

	.method( "Load", &Model::Load)
	.method( "Save", &Model::Save)
	;
}
