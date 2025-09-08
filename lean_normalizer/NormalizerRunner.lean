
import PhysicsExpr

open PhysicsExpr

/-- Main entry point -/
def main : IO Unit := do
  -- Simple test
  let exprs := [("rho", 1), ("z", 1), ("rho + z", 2)]
  let results := normalizeExpressionBatch exprs
  
  for info in results do
    IO.println s!"Expression: {info.expr}"
    IO.println s!"Signature: {info.signature}"
    IO.println s!"Depth: {info.depth}"
    IO.println s!"Normalized: {info.normalized}"
    IO.println s!"Valid: {info.isValid}"
    IO.println ""
