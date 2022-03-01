module Result = struct
  let (let*) = Result.bind
  let (let+) x f = Result.map f x
end

module Option = struct
  let (let*) = Option.bind
  let (let+) x f = Option.map f x
end

module Cmdliner_term = struct
  let (let+) x f = Cmdliner.Term.(const f $ x)
  let (and+) x y = Cmdliner.Term.(const (fun x y -> x, y) $ x $ y)
end
