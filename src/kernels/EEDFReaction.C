#include "EEDFReaction.h"

// MOOSE includes
#include "MooseVariable.h"

registerMooseObject("CraneApp", EEDFReaction);

template <>
InputParameters
validParams<EEDFReaction>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredCoupledVar("mean_energy", "The electron mean energy.");
  params.addRequiredCoupledVar("electrons", "The electron density.");
  params.addCoupledVar("target",
                       "The coupled target. If none, assumed to be background gas from BOLSIG+.");
  params.addRequiredParam<std::string>("reaction", "Stores the full reaction equation.");
  params.addRequiredParam<Real>(
      "coefficient", "The stoichiometric coefficient of this variable. (Gain or loss term.)");
  params.addParam<std::string>(
      "number",
      "",
      "The reaction number. Optional, just for material property naming purposes. If a single "
      "reaction has multiple different rate coefficients (frequently the case when multiple "
      "species are lumped together to simplify a reaction network), this will prevent the same "
      "material property from being declared multiple times.");
  return params;
}

EEDFReaction::EEDFReaction(const InputParameters & parameters)
  : Kernel(parameters),

    _reaction_coeff(getMaterialProperty<Real>("k" + getParam<std::string>("number") + "_" +
                                              getParam<std::string>("reaction"))),
    _d_k_d_actual_mean_en(getMaterialProperty<Real>("d_k" + getParam<std::string>("number") +
                                                    "_d_en_" + getParam<std::string>("reaction"))),
    _mean_en(coupledValue("mean_energy")),
    _em(coupledValue("electrons")),
    _mean_en_id(coupled("mean_energy")),
    _em_id(coupled("electrons")),
    _target(coupledValue("target")),
    _target_id(coupled("target")),
    _coefficient(getParam<Real>("coefficient"))
{
}

EEDFReaction::~EEDFReaction() {}

Real
EEDFReaction::computeQpResidual()
{
  return -_test[_i][_qp] * _em[_qp] * _target[_qp] * _reaction_coeff[_qp] * _coefficient;
}

Real
EEDFReaction::computeQpJacobian()
{
  if (_var.number() == _em_id)
  {
    Real d_k_d_em = _d_k_d_actual_mean_en[_qp] *
                    -(_mean_en[_qp] / (_em[_qp] * _em[_qp])) * _phi[_j][_qp];

    return -_test[_i][_qp] * _target[_qp] * _coefficient *
           (_reaction_coeff[_qp] * _phi[_j][_qp] + d_k_d_em * _em[_qp]);
  }
  else if (_var.number() == _target_id)
  {
    return -_test[_i][_qp] * _em[_qp] * _phi[_j][_qp] * _reaction_coeff[_qp] *
           _coefficient;
  }
  else
    return 0;
}

Real
EEDFReaction::computeQpOffDiagJacobian(unsigned int jvar)
{
  // Define multiplication factor (determining if product or reactant)

  Real d_k_d_mean_en = _d_k_d_actual_mean_en[_qp] * (1.0 / _em[_qp]) * _phi[_j][_qp];

  Real d_k_d_em = _d_k_d_actual_mean_en[_qp] *
                  -(_mean_en[_qp] / (_em[_qp] * _em[_qp])) * _phi[_j][_qp];

  if (jvar == _mean_en_id)
    return -_test[_i][_qp] * _em[_qp] * _target[_qp] * d_k_d_mean_en * _coefficient;

  //else if (jvar == _em_id)
  else if ((jvar == _em_id) && (_var.number() != _em_id))
    return -_test[_i][_qp] * _target[_qp] * _coefficient *
           (_reaction_coeff[_qp] * _phi[_j][_qp] + d_k_d_em * _em[_qp]);
  //else if (jvar == _target_id)
  else if ((jvar == _target_id) && (_var.number() != _target_id))
    return -_test[_i][_qp] * _em[_qp] * _phi[_j][_qp] * _reaction_coeff[_qp] *
           _coefficient;

  else
    return 0.0;
}
