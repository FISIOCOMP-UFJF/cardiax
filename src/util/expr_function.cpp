#include "util/expr_function.hpp"
#include "util/util.hpp"                 // show_error

// The ONLY place exprtk is included in the whole codebase. Keeping it in a
// single .cpp is what makes the build fast: exprtk is compiled once here
// instead of in every translation unit that includes poisson.hpp.
//
// exprtk has no `error` macro problem anymore (the util macro is show_error),
// so no push_macro/undef guard is needed. Enhanced features are disabled to
// trim exprtk's own (already heavy) compile time, since the RHS/exact-solution
// expressions only use +-*/ , ^ and the standard functions (sin, cos, ...).
#define exprtk_disable_enhanced_features
#include "thirdparty/exprtk.hpp"

struct ExprScalarFunction::Impl
{
  mutable double x = 0.0, y = 0.0, z = 0.0;
  exprtk::expression<double> expr;
};

ExprScalarFunction::ExprScalarFunction(const std::string & expr_str)
  : impl_(new Impl)
{
  exprtk::symbol_table<double> symbols;
  symbols.add_variable("x", impl_->x);
  symbols.add_variable("y", impl_->y);
  symbols.add_variable("z", impl_->z);
  symbols.add_constants();               // pi, e, epsilon

  impl_->expr.register_symbol_table(symbols);

  exprtk::parser<double> parser;
  if (!parser.compile(expr_str, impl_->expr))
    show_error(("ExprScalarFunction: cannot parse '" + expr_str + "'").c_str());
}

// Out-of-line destructor: required for the pimpl. It must be defined HERE,
// where Impl is complete, so unique_ptr<Impl> knows how to delete it. If it
// were defaulted in the header, the deleter would be instantiated where Impl
// is only forward-declared and fail to compile.
ExprScalarFunction::~ExprScalarFunction() = default;

double ExprScalarFunction::operator()(const arma::vec3 & p) const
{
  impl_->x = p(0);
  impl_->y = p(1);
  impl_->z = p(2);
  return impl_->expr.value();
}
