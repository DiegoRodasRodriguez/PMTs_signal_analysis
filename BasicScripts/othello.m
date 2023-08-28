function othello()

% othello.m
% A simple othello program
% Author: Subroto Gunawan
% Date: August 8, 2000

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% VARIABLE DECLARATIONS %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% drawing parameters
DrawResolution = 1/100;
N=2*pi;
Theta = 0:DrawResolution:N;

% define the board
BoardSize = 8;
FigSize = 400;
ButtonWidth = 100;
Board = zeros(BoardSize,BoardSize);
InitXLoc = [ BoardSize/2 BoardSize/2+1 BoardSize/2 BoardSize/2+1];
InitYLoc = [ BoardSize/2 BoardSize/2 BoardSize/2+1 BoardSize/2+1];

% default values
DefaultFirstTurn = 1; % 1 is white, 2 is black
finish=0;
whitePlayer = 'human';
blackPlayer = 'com';
XPos = [-1];
YPos = [-1];
XLegalPos = [];
YLegalPos = [];
BetweenCount = [];
NoMoveTurnCount = 0;
XDir = [-1 -1 0 +1 +1 +1 0 -1]; % counter clockwise directions starting from West
YDir = [0 +1 +1 +1 0 -1 -1 -1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% INITIALIZATION %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set up the main frame
fig1 = figure(1); 
clf;
set(fig1, 'Position', [120 120 FigSize+ButtonWidth FigSize],...
   'Name','Othello',...
   'NumberTitle','off',...
   'MenuBar','none');

% user control panel
left = 400;
bottom = 300;
width = 75;
height = 30;
ymod = 30;

statusText = uicontrol('Style','Text',...
   'Position',[(FigSize+ButtonWidth-300)/2 FigSize-30 300 25],...
   'BackgroundColor',[0.8 0.8 0.8],...
   'FontWeight','bold',...
   'String','');

infoText = uicontrol('Style','Text',...
   'Position',[(FigSize+ButtonWidth-300)/2 5 300 25],...
   'BackgroundColor',[0.8 0.8 0.8],...
   'String','Right Click to Quit');

whitePlayerText = uicontrol('Style','Text',...
   'Position',[left bottom width height],...
   'BackgroundColor',[0.8 0.8 0.8],...
   'HorizontalAlignment','left',...
   'String',[' White : ' whitePlayer]);
   
blackPlayerText = uicontrol('Style','Text',...
   'Position',[left bottom-ymod width height],...
   'BackgroundColor',[0.8 0.8 0.8],...
   'HorizontalAlignment','left',...
   'String',[' Black : ' blackPlayer]);

turnText = uicontrol('Style','Text',...
   'Position',[left bottom-4*ymod width height],...
   'BackgroundColor',[0.8 0.8 0.8],...
   'String',' White Turn ');
   
whiteText = uicontrol('Style','Text',...
   'Position',[left bottom-7*ymod width height],...
      'BackgroundColor',[0.8 0.8 0.8],...
      'String','White = 2');
   
blackText = uicontrol('Style','Text',...
   'Position',[left bottom-8*ymod width height],...
      'BackgroundColor',[0.8 0.8 0.8],...
      'String','Black = 2');

% creating the board grids
XLine = linspace(0,FigSize,BoardSize+1);
YLine = linspace(0,FigSize,BoardSize+1);

for i=1:(BoardSize+1)
   line([XLine(i) XLine(i)], [YLine(1) YLine(BoardSize+1)],'Color',[0 0 0]);
   line([XLine(1) XLine(BoardSize+1)], [YLine(i) YLine(i)],'Color',[0 0 0]);
end

axis([0,FigSize+ButtonWidth,0,FigSize]);
hold on; 
axis off;

% determine the radius of the piece
gap = (XLine(2)-XLine(1))/2;
radius = gap - (gap*0.25); % leave a 25% side margin

% save the center of the pieces for drawing purposes
XCenter = linspace((XLine(2)-XLine(1))/2, FigSize-(XLine(2)-XLine(1))/2, BoardSize);
YCenter = linspace((YLine(2)-YLine(1))/2, FigSize-(YLine(2)-YLine(1))/2, BoardSize);

BoardCenterX = repmat(XCenter,BoardSize,1);
BoardCenterY = repmat((FigSize - YCenter)',1,BoardSize);

% set up the transparent pieces and draw the first four pieces
Board(InitYLoc(1),InitXLoc(1)) = 1;
Board(InitYLoc(2),InitXLoc(2)) = 2;
Board(InitYLoc(3),InitXLoc(3)) = 2;
Board(InitYLoc(4),InitXLoc(4)) = 1;

for i=1:BoardSize
   for j=1:BoardSize
      X = BoardCenterX(j,i) + radius*cos(Theta);
      Y = BoardCenterY(j,i) + radius*sin(Theta);
      switch Board(j,i)
         case 1, BoardColor(j,i) = fill(X,Y,'w');
         case 2, BoardColor(j,i) = fill(X,Y,'k');
         otherwise
         	BoardColor(j,i) = fill(X,Y,[0.8 0.8 0.8]);
            set(BoardColor(j,i),'EdgeColor',[0.8 0.8 0.8]);
      end
   end
end

%%%%%%%%%%%%%%%%%%%%%%
%%%% MAIN PROGRAM %%%%
%%%%%%%%%%%%%%%%%%%%%%

turn = DefaultFirstTurn;

while ~finish
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % check and see if there is any move available
      
      canMove = 0;
      XPosLegal = [];
      YPosLegal = [];
      BetweenCount = [];

      for i=1:BoardSize
         for j=1:BoardSize
            
            if Board(j,i)==0
                              
               for dir=1:8
                  
                  done=0;
                  count=2;
                  between=[];
                                    
                  while(~done)
                     markerX = i + count*XDir(dir);
                     markerY = j + count*YDir(dir);
                     if (markerX<1) | (markerX>BoardSize) | (markerY<1) | (markerY>BoardSize)
                        done=1;
                     elseif (Board(markerY,markerX)==turn)
                        done=1;
                        between = [between Board(markerY-YDir(dir),markerX-XDir(dir))];
                        
                        if (Board(markerY,markerX)==turn) & (between==(mod(turn,2)+1))
                           canMove=1;
                           XPosLegal = [XPosLegal i];
                           YPosLegal = [YPosLegal j];
                           BetweenCount = [BetweenCount length(between)];
                        end
                     else
                        count = count+1;
                        between = [between Board(markerY-YDir(dir),markerX-XDir(dir))];
                     end
                     
                  end % end of while not done
                  
               end % end of for dir
               
            end % end of if board==0
            
         end
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % A piece of code for debugging purposes
      
      %BoardLegal = zeros(BoardSize,BoardSize);
      %for i=1:length(XPosLegal)
      %   BoardLegal(YPosLegal(i),XPosLegal(i))=1;
      %end
      %BoardLegal
      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % if there is a move available and white turns, get the input from player
      
if canMove
         
   if turn==1; % if it's the player's turn
   
      done=0;
      while (~done)
         
         [Xnew Ynew Button] = ginput(1);
         if Button==1
            
            if (XLine(1)<=Xnew) & (Xnew<=XLine(BoardSize+1)) & ...
                  (YLine(1)<=Ynew) & (Ynew<=YLine(BoardSize+1) ) 
               for i=1:BoardSize
                  if (XLine(i)<Xnew) & (Xnew<XLine(i+1))
                     XPos = i;
                  end
                  if (YLine(i)<Ynew) & (Ynew<YLine(i+1))
                     YPos = BoardSize+1-i;
                  end
               end
            end
         
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % based on previous result, check and see if the move is legal
         % and if it's legal, draw the piece
         
         legal = 0;
         for i=1:length(XPosLegal)
            if (XPos == XPosLegal(i)) & (YPos == YPosLegal(i))
               legal = 1;
            end
         end
                  
         if legal
            NoMoveTurnCount = 0;
            Board(YPos,XPos) = turn;
            set(BoardColor(YPos,XPos),'FaceColor','w','EdgeColor',[0 0 0]);
            set(statusText,'String','');
            done=1;
            drawnow;
         else
            set(statusText,'String','***** Invalid Move *****');
            drawnow;
         end
               
         else                % if right click
            done=1;
            finish=1;
         end
      end % end of while not done
   
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % then check into 8 direction, reverse the pieces and draw
   
   if ~finish
      
   for dir=1:8
      done=0;
      count=2;
      between=[];
      
      while(~done)
         markerX = XPos + count*XDir(dir);
         markerY = YPos + count*YDir(dir);
         if (markerX<1) | (markerX>BoardSize) | (markerY<1) | (markerY>BoardSize)
            done=1;
         elseif (Board(markerY,markerX)==turn)
            done=1;
            between = [between Board(markerY-YDir(dir),markerX-XDir(dir))];
            if (Board(markerY,markerX)==turn) & (between==(mod(turn,2)+1))
               for i=1:(count-1)
                  Board(YPos+YDir(dir)*i,XPos+XDir(dir)*i) = turn;
                  set(BoardColor(YPos+YDir(dir)*i,XPos+XDir(dir)*i),...
                     'FaceColor','w','EdgeColor',[0 0 0]);
                  drawnow;
               end
            end
         else
            count = count+1;
            between = [between Board(markerY-YDir(dir),markerX-XDir(dir))];
         end
      end
   end
  
   end % end of if ~finish
   
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % if it is black turn, activate the AI
   
else
   % search the highest number of pieces reversed
   [dummy Pos] = max(BetweenCount);
   XPos = XPosLegal(Pos);
   YPos = YPosLegal(Pos);
   Board(YPos,XPos) = turn;
   
      
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % draw the new piece
   NoMoveTurnCount = 0;
   set(statusText,'String','');
   set(BoardColor(YPos,XPos),'FaceColor','k','EdgeColor',[0 0 0]);
   drawnow;
   
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % then check into 8 direction, reverse the pieces and draw
   
   for dir=1:8
      done=0;
      count=2;
      between=[];
      while(~done)
         markerX = XPos + count*XDir(dir);
         markerY = YPos + count*YDir(dir);
         if (markerX<1) | (markerX>BoardSize) | (markerY<1) | (markerY>BoardSize)
            done=1;
         elseif (Board(markerY,markerX)==turn)
            done=1;
            between = [between Board(markerY-YDir(dir),markerX-XDir(dir))];
            if (Board(markerY,markerX)==turn) & (between==(mod(turn,2)+1))
               for i=1:(count-1)
                  Board(YPos+YDir(dir)*i,XPos+XDir(dir)*i) = turn;
                  set(BoardColor(YPos+YDir(dir)*i,XPos+XDir(dir)*i),...
                     'FaceColor','k','EdgeColor',[0 0 0]);
                  drawnow;
               end
            end
         else
            count = count+1;
            between = [between Board(markerY-YDir(dir),markerX-XDir(dir))];
         end
      end
   end
   
end % end of separating human player and com player
   
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % if there is no legal move, display the status
   
else
   % this situation is only allowed twice, otherwise it will repeat forever
   if NoMoveTurnCount<3
      set(statusText,'String','***** No Legal Moves, Change Turns *****');
      drawnow;
      NoMoveTurnCount = NoMoveTurnCount + 1;
   else
      set(statusText,'String','***** No Legal Move For Both Players, End of Game *****');
      finish=1;
   end
end
   
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % once all is done, change turns
   
   turn = mod(turn,2)+1;
   if ~finish
   	if turn==1
      	set(turnText,'String','White Turn');
   	else
      	set(turnText,'String','Black Turn');
      end
   end
      
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % count the number of pieces and see if finished
   
   White = sum(sum(Board == ones(BoardSize,BoardSize)));
   Black = sum(sum(Board == repmat([2],BoardSize,BoardSize)));
   set(whiteText,'String',['White = ' num2str(White)]);
   set(blackText,'String',['Black = ' num2str(Black)]);
   
   if (White+Black) == BoardSize^2
      finish = 1;
      if (White > Black), winner = 'White';
      elseif (White < Black), winner = 'Black';
      else winner = 'Stalemate, no one';
      end
      
      set(statusText,'String',['***** ' winner ' wins! *****']);
   end
      
end % end of while - end of game