module Result = struct
  let (let*) = Result.bind
  let (let+) x f = Result.map f x
end

module Option = struct
  let (let*) = Option.bind
  let (let+) x f = Option.map f x
end

