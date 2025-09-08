import Lake
open Lake DSL

package «physics_normalizer» {
  -- add package configuration here
}

require mathlib from git
  "https://github.com/leanprover-community/mathlib4.git"

@[default_target]
lean_lib «PhysicsExpr» {
  -- add library configuration here
}

lean_lib «FoliationValidator» {
  -- Foliation constraint validator
}

@[default_target]
lean_exe «normalizer» {
  root := `NormalizerRunner
}

lean_exe «foliation_search» {
  root := `FoliationSearch
}

lean_exe «foliation_search_depth4» {
  root := `FoliationSearchDepth4
}