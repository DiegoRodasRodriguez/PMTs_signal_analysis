function [Out1,Out2,Out3,Out4,Out5,Out6,Out7,Out8,Out9,Out10,Out11,Out12]=cut(I,In1,In2,In3,In4,In5,In6,In7,In8,In9,In10,In11,In12);

%I=eval(expr);

%length(I)

Out1=In1(I);
if nargin>2,Out2=In2(I); end
if nargin>3,Out3=In3(I); end
if nargin>4,Out4=In4(I); end
if nargin>5,Out5=In5(I); end
if nargin>6,Out6=In6(I); end
if nargin>7,Out7=In7(I); end
if nargin>8,Out8=In8(I); end
if nargin>9,Out9=In9(I); end
if nargin>10,Out10=In10(I); end
if nargin>11,Out11=In11(I); end
if nargin>12,Out12=In12(I); end

return

