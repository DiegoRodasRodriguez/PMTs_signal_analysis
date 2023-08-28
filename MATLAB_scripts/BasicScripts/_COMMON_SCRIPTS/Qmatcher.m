function NewQm = Qmatcher(Qleft, Qright, e, e_match, lambda)

WeightLeft  = (1 - tanh(lambda * (e - e_match))) / 2;
WeightRight = (1 + tanh(lambda * (e - e_match))) / 2;

NewQright = Qright .* WeightRight;
NewQleft  = Qleft  .* WeightLeft;

NewQm     = NewQleft + NewQright;

end

