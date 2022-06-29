classdef MouseTrackerKF < MouseTracker
    properties
        default_thresh = 0;
        used_thresh = 0;
        p_mouse = .0007; %the proportion of pixels that are 
        kf = struct('nstates', {}, 'nob', {}, 's', {});
        %mm_conv = .862; %mm/px linear, old res for data pre 2015
        mm_conv = .7143; %mm/px linear, new resolution (14px/10mm) 
        fullPaths; %the full detected paths, to be saved if the paths are skeletonized
        exploredProp; %the proportions of the trails explored at each timepoint
        exploredLen; %the lengths of the trails explored at each timepoint (in px)
        trailShadow = [];
        trailShadow_size = 20;
        pathVerticesDefined = 0;
        pathVertices;
        colorBG = [];
    end
    
    methods
        function this = MouseTrackerKF(varargin)
            % MouseTrackerKF(varargin) - This is an object that loads a video and then detects a mouse in that video
            % letting you get the mouse position in various parts of the video.  It does this on an request basis, then 
            % saves the tracking results to minimize memory load. Also, this code is set up to track the position of
            % multiple portions of a binary thresholded image. This is meant to be different parts of the mouse when
            % viewed from below (feet, nose, and tail). Thus, the statistics compiled assume that the whole set of
            % "blobs" represents a single mouse.
            % 
            % Args: 1 - A filename or folder of the movie to load.  If there isn't one given, then a dialog is
            % presented.
            % 2 - frame range - a vector of frame numbers to consider as the whole movie, example 1:100
            % 3 - time range - a beginning and end time to consider as the whole movie, in sec, example [60 120].
            
            % This is all to just make sure that the file is opened correctly, and that it will take a filename,
            % directory or nothing.
            
            % All of the proper initialization is done in the MouseTracker object this is inherited from
%             disp('Initializing MouseTrackerKF');
%             if ~exist('filename', 'var') || isempty(filename)
%                 disp('Filename not given, exiting');
%                 this = [];
%                 return;
%             end
            this = this@MouseTracker(varargin{:});
            %this.blobID = NaN*zeros(size(this.areas)); %implementing a unique ID assigment for each blob to facilitate assignment
            
        end %function MouseTrackerKF
        % ------------------------------------------------------------------------------------------------------
        function initKalmanFilter(this)
        % function initKalmanFilter(this)
        %
        % Initialization of Kalman filter parameters
            this.kf = struct('nstates', {}, 'nob', {}, 's', {});
            zs = 4;
            this.kf(1).nstates = zs;
            this.kf(1).nob = 2;
            this.kf(1).s.A = eye(zs);
            this.kf(1).s.A = [1 0 1 0; 0 1 0 1; 0 0 1 0; 0 0 0 1]; 
            this.kf(1).s.z = nan*zeros(zs,1); 
            this.kf(1).s.x = nan*zeros(zs,1); 
            this.kf(1).s.H = eye(zs);
            this.kf(1).s.R = eye(zs); %measurement error covariance
            %this.kf.s.P = cov(np');
            this.kf(1).s.P = eye(zs);
            this.kf(1).s.u = zeros(zs,1);
            this.kf(1).s.B = eye(zs);
            this.kf(1).s.Q = eye(zs);

        end
        % ------------------------------------------------------------------------------------------------------
        function [wantedCOM, nose] = mousePosition(this, time_range)
        % function [wantedCOM, nose] = mousePosition(this, time_range)    
         % Returns the position of the mouse center of mass (body) and 
         % the nose for the specified time range
            wantedFrames = this.timesToFrames(time_range); 
            wantedCOM = this.bodyCOM(wantedFrames,:);
            if isnan(wantedCOM)
                nanFrames = wantedFrames(isnan(wantedCOM(:,1)));
                this = this.findMouse(nanFrames);
            end 
            % get what they wanted
            wantedCOM = this.bodyCOM(wantedFrames,:);
            nose = this.nosePos(wantedFrames, :);
        end
        
        % ------------------------------------------------------------------------------------------------------
        function vel = mouseVelocity(this, time_range)
            % function vel = mouseVelocity(this, time_range)
            % Returns the velocity of the mouse body for the specified frames
            frames = this.timesToFrames(time_range);
            vel = this.vel(frames);
            if isnan(vel)
                this.mousePosition(time_range);
            end
            vel = this.computeVelocity(frames);
        end
        
        % ----------------------------------------------------------------------------------------------------
        function tailPos = tailPosition(this, frames)
        % function tailPos = tailPosition(this, frames)
        %
        % We aren't going to want the tail position all that much, so assemble it each time
            tv = logical(this.tailVisible(frames));
            tb = this.tailblob(frames);
            tailPos = NaN*zeros(length(frames), 2);
            areas = this.areas(frames,:);
            for ii = 1:length(areas)
                if tv(ii)
                    tailPos(ii,:) = areas(ii,tb(ii)).Centroid;
                end
            end
        end
        
        % ------------------------------------------------------------------------------------------------------
        function ah = plotVelocity(this, frames, varargin)
        % function ah = plotVelocity(this, frames, varargin)
        %
        % Plots the position with the velocity as a color code. Varargin provides the following
        % optional inputs:
        % 1) 'nose' or 'body' for different velocities.  Nose is default.
        % 2) Gaussian filter length in samples - the velocity can be noisy, which can interfere with
        %    its visualization.  Filtering to reduce extremes.  Default is 1, the std in samples.
            mm_conv = .862; %mm/px linear
            
            %input parsing/checking
            if ~exist('frames', 'var') || isempty(frames); frames = 1:this.nFrames; end
            if length(varargin) > 0
                whichVel = varargin{1};
            else
                whichVel = 'nose';
            end
            if length(varargin) > 1
                filt_std = varargin{2};
            else
                filt_std = 1;
            end
            
            if strcmp(whichVel, 'nose')
                vel_vect = this.noseVel(frames);
            else
                vel_vect = this.bodyVel(frames,:);
                vel_vect = sqrt(sum(vel_vect.^2,2));
            end
            vel_vect = gaussianFilter(vel_vect, filt_std, 'conv');
            vel_vect = mm_conv * this.frameRate * vel_vect;  %convert to mm/sec
            [cm cinds] = getIndexedColors('jet', vel_vect, 0, [0 300]);
            
            % unfortunately, I haven't been able to figure out an easier way to do colormapping of points plotted over
            % a B&W image, without the points and the image sharing the same colormap.  So, here I'm just plotting every
            % point separately, with a color indicative of the velocity of the animal.
            figure;
            imshow(this.plotPathsOnBg()); hold on;
            xlim([0 this.width]); %fit the axes to the image
            ylim([0 this.height]);
            for ii=1:length(frames)
                fi = frames(ii);
                %ci = find(sorted_vel == this.vel(fi), 1, 'first');
                if ~isnan(vel_vect(ii))
                    line('Xdata', this.nosePos(fi,1), 'Ydata', this.nosePos(fi,2), 'Marker', '.', 'MarkerSize', 12, 'Color', cm(cinds(ii),:));
                end
            end 
            ah = gca;
            %Make another figure entirely to get a color scale
            fh2 = figure;
            %ah2 = axes('position', [.1 .1 .7 .8], 'visible', 'off');
            colormap(cm);
            %axes(ah2);
            %pcolor([vel_vect(frames,1), vel_vect(frames,1)]);
            pcolor([0:300; 0:300]);
            %line('parent', ah2, 'ydata', sorted_dir, 'xdata', 1:length(sorted_dir)); 
            colorbar;
        end
        % ------------------------------------------------------------------------------------------------------
        function plotVelocityTimes(this, time_range)
         % Version using a range of times, rather than frames
            frames = this.timesToFrames(time_range);
            this.plotVelocity(frames);
        end
        
        % ------------------------------------------------------------------------------------------------------
        function plotPosition(this, frames, varargin)
            % function plotPosition(this, frames)
            %
            % This plots the detected mouse positions in red over the background frame
            c = [.3, .6, 1]; %plotblue
            marker = '.';
            ah = [];
            overwrite = 1;
            if length(varargin) >=1
                ah = varargin{1};
                newFig = 0;
            end
            if length(varargin) >= 2
                overwrite = varargin{2};
            end
            if length(varargin) >=3
                c = varargin{3};
            end
            if length(varargin) >=4
                marker = varargin{4};
            end
            if isempty(frames); frames = 1:this.nFrames; end
            if isempty(ah)
                fh = figure; 
                ah = axes('position', [.1, .1, .7 .8]); hold on;
                newFig = 1;
            end
            if (newFig || overwrite) 
                pathIm = this.plotPathsOnBg();
                imshow(pathIm); hold on;
            end
            xlim([0 this.width]);
            ylim([0 this.height]);
            for ii=1:length(frames)
                fi = frames(ii);
                %line('Parent', ah, 'Xdata', this.bodyCOM(fi,1), 'Ydata', this.bodyCOM(fi,2), ...
                %       'Marker', '.','MarkerSize', 8, 'Color', c);
                if ~isnan(this.noseblob(fi))
                    % So, this way we can't plot a solid line.
                    line('Parent', ah, 'Xdata', this.nosePos(fi,1), 'Ydata', this.nosePos(fi,2), ...
                        'Marker', '.','MarkerSize', 8, 'Color', c);
                end
            end
            
        end
        % ------------------------------------------------------------------------------------------------------
        function plotFilterPosition(this,frames)
        % function plotFilterPosition(this,frames)
            plotblue = [.3, .6, 1]; % a nice blue color
            if ~exist('frames', 'var') || isempty(frames); frames = 1:this.nFrames; end
            fh = figure;
            ah = axes('position', [.1, .1, .7 .8]);
            bgIm = this.plotPathsOnBg();
            imshow(bgIm); hold on;
            xlim([0 this.width]);
            ylim([0 this.height]);
            filt_x = [this.kf.s.x]; %the kalman filter positions
            filt_x = filt_x(1:2,:)';
            for ii=1:length(frames)
                fi = frames(ii);
                if ~isnan(this.noseblob(fi))
                    line('Parent', ah, 'Xdata', this.areas(fi,this.noseblob(fi)).Centroid(1), 'Ydata', ...
                        this.areas(fi,this.noseblob(fi)).Centroid(2), 'Marker', '.', 'MarkerSize', 8, 'Color', plotblue);
                    line('Parent', ah, 'Xdata', filt_x(fi,1), 'Ydata', filt_x(fi,2), 'Marker', '.', 'MarkerSize', 8, 'Color', [1 0 0]);
                end
            end 
            title(this.videoFN);
        end
        % ------------------------------------------------------------------------------------------------------
        function plotNosePosition(this, frames)
            % function plotDirection(this, frames)
            plotblue = [.3, .6, 1]; % a nice blue color
            if ~exist('frames', 'var') || isempty(frames); frames = 1:this.nFrames; end
            fh = figure;
            ah = axes('position', [.1, .1, .7 .8]);
            bgIm = this.plotPathsOnBg();
            %bgIm = this.plotPaths();
            imshow(bgIm); hold on;
            xlim([0 this.width]);
            ylim([0 this.height]);
            for ii=1:length(frames)
                fi = frames(ii);
                if ~isnan(this.noseblob(fi))
                    np = this.nosePos(fi,:);
                    %line('Parent', ah, 'Xdata', this.areas(fi,this.noseblob(fi)).Centroid(1), 'Ydata', ...
                    %    this.areas(fi,this.noseblob(fi)).Centroid(2), 'Marker', '.', 'MarkerSize', 8, 'Color', plotblue);
                    line('Parent', ah, 'Xdata', this.nosePos(fi,1), 'Ydata', ...
                        this.nosePos(fi,2), 'Marker', '.', 'MarkerSize', 10, 'Color', plotblue);
%                     if(mod(ii, 50) == 0)
%                        text(np(1)+5, np(2)+5, num2str(fi));
%                     end
                    %line('Parent', ah, 'Xdata', this.nosePos(fi,1), 'Ydata', ...
                    %    this.nosePos(fi,2), 'Marker', '.', 'MarkerSize', 10, 'Color', 'w');
                end
            end 
            title(this.videoFN);
        end
        % ------------------------------------------------------------------------------------------    
        function [nb_stat, nose_bright] = getNoseBrightness(this, frames, operation)
            % operation is a function pointer to operation to be performed
            % on the set of nose pixels 
            nSegs = ceil(length(frames)/this.framesPerSeg);
            movieDone = 0; 
            dbg = 0;
            nose_bright = NaN*zeros(length(frames), this.MAX_SIZE_THRESH);
            nb_stat = NaN*zeros(length(frames),1);
            for jj = 1:nSegs %divide the computation up into segments so that there is never too much in memory at once
                first = (jj-1)*this.framesPerSeg + 1; %first frame of segment
                if (jj == nSegs) % the last frame
                    last = length(frames);
                else
                    last = (jj)*this.framesPerSeg;
                end
                disp(sprintf('Finding brightnesses for frames %d - %d', first, last));
                segFrames = frames(first:last);
                
                frameArray = this.readMovieSection(segFrames,'diff');
                for ii = 1:(last-first)
                    fi = frames(first+ii-1);
                    if ~isnan(this.noseblob(fi))
                        ni = this.areas(fi, this.noseblob(fi)).PixelIdxList;
                        frame = frameArray(:,:,ii);
                        nose_px = frame(ni);
                        nose_bright(first+ii-1, 1:length(nose_px)) = nose_px; 
                        nb_stat(first+ii-1) = operation(nose_px(:));
                    end
                end
                        
            end
        end
        % -------------------------------------------------------------------------------------------
        function plotDirection(this, frames)
            % function plotDirection(this, frames)
            plotblue = [.3, .6, 1]; % a nice blue color
            if ~exist('frames', 'var') || isempty(frames); frames = 1:this.nFrames; end
            fh = figure;
            ah = axes('position', [.1, .1, .7 .8]);
            sorted_dir = sort(this.direction); %sorted vector for colormapping
            [cm cinds] = getIndexedColors('jet', this.direction(frames), 1);
            bgIm = this.plotPathsOnBg();
            imshow(bgIm); hold on;
            xlim([0 this.width]);
            ylim([0 this.height]);
            % unfortunately, I haven't been able to figure out an easier way to do colormapping of points plotted over
            % a B&W image, without the points and the image sharing the same colormap.  So, here I'm just plotting every
            % point separately, with a color indicative of the velocity of the animal.
            for ii=1:length(frames)
                fi = frames(ii);
                ci = find(sorted_dir == this.direction(fi), 1, 'first');
                if ~isnan(this.noseblob(fi))
                    line('Xdata', this.nosePos(fi,1), 'Ydata', ...
                        this.nosePos(fi,2), 'Marker', '.', 'MarkerSize', 8, 'Color', cm(cinds(ii),:));
                end
                if ~isnan(ci)
                    %line('Xdata', this.COM(fi,1,1), 'Ydata', this.COM(fi,2,1), 'Marker', '.', 'MarkerSize', 10, 'Color', cm(ci,:));
                    %line('Xdata', this.COM(fi,1,1), 'Ydata', this.COM(fi,2,1), 'Marker', '.', 'MarkerSize', 10, 'Color', 'r');
                    %line('parent', ah2, 'Xdata', this.COM(ii,1), 'Ydata', this.COM(ii,2), 'Marker', '.', 'MarkerSize', 10, 'Color', cm(ci,:));
                end
            end 
            title(this.videoFN);
            % Make another graph with the color scale
            fh2 = figure;
            ah2 = axes('position', [.1 .1 .7 .8], 'visible', 'off');
            colormap(cm);
            axes(ah2);
            pcolor(repmat(sorted_dir(:), [1 3]));
%             line('parent', ah2, 'ydata', sorted_dir, 'xdata', 1:length(sorted_dir)); 
            colorbar;
        end
        
        % ---------------------------------------------------------------------------
        function plotFollowing(this, frames, dist_thresh, textflag)
            % function plotFollowing(this, frames, dist_thresh, textflag)
            % plot the following segments based on the following information 
            if isempty(frames)
                frames = 1:this.nFrames;
            end
            [dists, fframes] = this.distanceOnTrail(frames,1,dist_thresh);
            this.plotPosition(frames); hold on;
            for i=1:size(fframes,1)
                range = fframes(i,1):fframes(i,2);
                np = this.nosePos(range, :);
                plot(np(:,1), np(:,2), '.m', 'MarkerSize',10);
                if textflag
                    text(mean(np(:,1)), mean(np(:,2)), num2str(dists(i)), 'Color', 'w', 'FontSize', 12);
                end
            end
        end
        
        % ---------------------------------------------------------------------------
        function ah = plotFollowingSide(this, frames, dist_thresh, textflag)
            % function ah = plotFollowingSide(this, frames, dist_thresh, textflag)
            % plot the following segments based on the following information
            if isempty(frames)
                frames = 1:this.nFrames;
            end
            [dists, fframes] = this.distanceOnTrail(frames,1,dist_thresh);
            this.plotNosePosition(frames); hold on;
            for i=1:size(fframes,1)
                range = fframes(i,1):fframes(i,2);
                np = this.nosePos(range, :);
                %incorporate the side of the trail coloring
                signed_dists = this.orthogonalDistFromTrail(range,1);
                negd = signed_dists < 0;
                posd = signed_dists > 0;
                plot(np(posd,1), np(posd,2), '.m', 'MarkerSize',10);
                plot(np(negd,1), np(negd,2), '.y', 'MarkerSize',10);
                if textflag
%                     for jj = 1:length(range)
%                         text(np(jj,1), np(jj,2), num2str(signed_dists(jj)), 'Color', 'w');
%                     end
                end
            end
            ah = gca;
        end
        % ---------------------------------------------------------------------------
        function plotOrientation(this, frames)
            % function plotOrientation(this, frames)
            if ~exist('frames', 'var') || isempty(frames); frames = 1:this.nFrames; end
            %orient = this.headingFromBodyToNose(frames);
            orient = this.headingFromMotion_TrailRelative(frames,1);
            cm_vals = linspace(-pi,pi,128);
            [cm, ~, cvals] = getIndexedColors('colormapc', cm_vals, 0);
            cinds = zeros(size(orient))*NaN;
            for ii=1:length(orient)
                dif = abs(cvals - orient(ii));
                [~, di] = min(dif);
                if ~isnan(di)
                    cinds(ii) = di;
                end
            end
            figure;
            imshow(this.plotPathsOnBg()); hold on;
            xlim([0 this.width]); %fit the axes to the image
            ylim([0 this.height]);
            for ii=1:length(frames)
                fi = frames(ii);
                %ci = find(sorted_orient == this.bodyOrient(fi), 1, 'first');
                if ~isnan(orient(ii))
                    line('Xdata', this.nosePos(fi,1), 'Ydata', this.nosePos(fi,2), 'Marker', '.', 'MarkerSize', 10, 'Color', cm(cinds(ii),:));
                end
            end 
            title(this.videoFN);
            % Make another graph with the color scale
            fh2 = figure;
            ah2 = axes('position', [.1 .1 .7 .8]);
            makeCircColorbar(0,cm,cvals,0);
            set(gca, 'XDir', 'reverse')
        end
        % ---------------------------------------------------------------------------------------------------
        function [exploredProp, exploredLen] = plotFollowingTimecourse(this, threshDist, varargin)
            % function plotFollowingTimecourse(this)
            %
            % This function returns the proporation of the trail pixels that the mouse came within a given
            % distance of, for each video frame.  It's nice because it does show the bouts of following in a 
            % way that even penalizes the animal for stopping along the trail. It's a  clean way of measuring 
            % the completeness of his trail exploration.
            % Algorthmically, this is a similar problem to finding the following segments except reversed.  As
            % implemented it's SLOW.  Thought for speeding it up is to do it on the full video once, to get the
            % frames where the animal is within distance of ANY pixel, then do the incremental on only those frames.
            %
            % exploredLen returns a length in mm using mm_conv
            pb=1;
            if (sum(sum(this.exploredProp)) > 0)
                calculated = 1;
            else 
                calculated = 0;
            end
            if nargin > 2
                ah = varargin{1};
                if isempty(ah)
                    pb = 0;
                end
            else
                % plotting part
                figure;
                ah = axes; hold on;
            end
            if ~calculated
                nPaths = length(this.paths);
                exploredProp = zeros(this.nFrames, nPaths);
                exploredLen = zeros(this.nFrames, nPaths);
                this.makePathsSkel();
                for ii = 1:this.nFrames
                    if ~mod(ii, 100)
                        disp('Completed 100 frames');
                    end
                    np = this.nosePos(1:ii,:);
                    nn = ~isnan(np(:,1));
                    np = np(nn,:);
                    for trailNum = 1:nPaths
                        trailPos = this.paths(trailNum).PixelList;
                        
                        distm = ipdm(single(np), single(trailPos));
                        [trailDist, mini] = nanmin(distm, [], 1);
                        explored = find(trailDist <= threshDist);
                        npx = length(this.paths(trailNum).PixelIdxList);
                        exploredProp(ii,trailNum) = length(explored)/npx;
                        exploredLen(ii,trailNum) = length(explored);
                    end
                end
                this.makePathsFull();
                this.exploredProp = exploredProp;
                this.exploredLen = exploredLen * this.mm_conv;
            else
                exploredProp = this.exploredProp;
                exploredLen = this.exploredLen;
            end
            
            % plotting part
            if pb
                plot(this.times/1000, this.exploredLen(:,1), 'g','LineWidth', 2); hold on;
                plot(this.times/1000, this.exploredLen(:,2), 'r','LineWidth', 2); 
                xlabel('Time (sec)');
                ylabel('Length of the trail explored');
            end
            
        end
        
       
        % ------------------------------------------------------------------------------------------------------
        function writeMovie(this, filename, movieType, frames, dispCrop, speed)
            % function writeMovie(this, filename, movieType, frames, dispCrop)
            % Writes a movie of the tracked mouse to disk at the given location
            if isempty(frames)
                frames = 1:this.nFrames;
            end
            vidWriter = VideoWriter(filename, 'MPEG-4');
            
            vidWriter.FrameRate = round(speed*this.frameRate);
            %vidWrite.Quality = 100;
            open(vidWriter);
            figure;
            for ii=frames
                this.showFrame(ii, movieType, dispCrop);
                currFrame = getframe; % now to the business of writing the movie
                writeVideo(vidWriter,currFrame);
                hold off;
            end
            close(vidWriter);
        end
        
        % -------------------------------------------------------------------------------------------------
        function showMovie(this, movieType, frames, varargin)
            % function showMovie(this, useBinMovie, frames, dispCrop[OPTIONAL], overlay[OPTIONAL])
            if nargin >=4
                dispCrop = varargin{1};
            else
                dispCrop = [];
            end
            if nargin >=5
                overlayMov = varargin{2};
            end
            fh = figure;
            this.exitMovie = 0;
            set(fh, 'WindowKeyPressFcn', @this.exitMovieLoop);
            if ~exist('frames','var') || isempty(frames)
                frames = 1:this.nFrames;
            end
            for ii=1:length(frames)
                fi = frames(ii);
                if this.exitMovie
                    break
                end
                if nargin < 5
                    this.showFrame(fi, movieType, dispCrop);
                else
                    this.showFrame(fi, movieType, dispCrop, overlayMov(:,:,ii));
                end
                hold off;
                pause(1/this.frameRate/10); %faster than the natural framerate due to impatience.
            end
            this.exitMovie = 0;
            set(fh, 'WindowKeyPressFcn', '');
        end
        
        % -------------------------------------------------------------------------------------------------
        function exitMovieLoop(this, src, event)
            % Function to set a flag internally to exit a movie that is being displayed.
            % It is not used for any other purpose
            if strcmp(event.Key, 'q') || strcmp(event.Key, 'escape')
                this.exitMovie = 1;
            end
        end
        
        function mov = returnFrames(this, frames, movieType)
            % mov = returnFrames(this, frames, binary)
            % returns the frames of the movie specified, binary or greyscale
            mov = this.readMovieSection(frames, movieType);
            
        end
        
        % Setting functions for manually altering the tracking result
        function setTailPosition(this, frame, pos)
            frame = frame(1);
            centers = [this.areas(frame, :).Centroid];
            centers = reshape(centers, 2, [])';
            dist = ipdm(centers, pos);
            if nanmin(dist) < 10
                [min_dist, disti] = nanmin(dist);
                this.tailblob(frame) = disti;
            end
        end
        
        function setNosePosition(this, frame, pos)
            frame = frame(1);
            centers = [this.areas(frame, :).Centroid];
            centers = reshape(centers, 2, [])';
            dist = ipdm(centers, pos);
            if nanmin(dist) < 15
                [min_dist, disti] = nanmin(dist);
                this.noseblob(frame) = disti;
                this.nosePos(frame,:) = this.areas(frame,disti).Centroid;
            end
        end
        
    %end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRIVATE METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %methods %(Access = private)
        
        function this = findMouse(this, frames)
            % function this = findMouse(this, frames)
            
            nSegs = ceil(length(frames)/this.framesPerSeg);
            movieDone = 0; 
            
            for jj = 1:nSegs %divide the computation up into segments so that there is never too much in memory at once
                first = (jj-1)*this.framesPerSeg + 1; %first frame of segment
                if (jj == nSegs) % the last frame
                    last = length(frames);
                else
                    last = (jj)*this.framesPerSeg;
                end
                segFrames = frames(first:last);
                disp(['Finding mouse in segment ' num2str(jj)]);
                
                if ~this.logFile
                    [frameArray, fcLum] = this.readMovieSection(segFrames,'bin', this.default_thresh);
                else
                    frameArray = this.makeMovieFromLog(segFrames); 
                    fcLum = zeros(length(segFrames),1);
                end
                segFrames = segFrames(1:size(frameArray, 3)); % look at the return to verify it returned the correct size array
                this.fcLum(segFrames) = fcLum; % save the frame counter area value
                newFrameCount(jj) = this.nFrames;
                trimmed(jj) = 0; 
                if size(frameArray,3) < length(segFrames) %there was a problem reading some frames, so we need to adjust
                    disp('In MouseTracker.findMouse, frameArray has too few frames, revising expectations');
                    fdiff = this.nFrames - segFrames(end); %now, consider the entire length of the frames
                    newFrameCount(jj) = this.nFrames - fdiff; 
                    trimmed(jj) = 1;
                    %movieDone = 1; %don't read movie anymore
                end
                detectMouseInSection(this, segFrames, frameArray);
                clear frameArray;
            end
            if trimmed %means the last set of frames came out shorter than expected from movie info, trimming things
                this.nFrames = min(newFrameCount);
                frames = frames(1:this.nFrames); %just trim them off the end - this error only happens for last seg
                this.trimFields(); %trim the object property arrays
            end
            this.computeVelocity(frames);
            this.refineTracking(frames);    
        end
        
        % -----------------------------------------------------------------------------------------------
        function detectMouseInSection(this, segFrames, frameArray)
            
            dbg = 0;
            cb = 0;
            
            for ii = 1:length(segFrames)  %have to work 1 frame at a time, unfortunately
                lm = bwlabel(frameArray(:,:,ii)); %label matrix
                fi = segFrames(ii);
                temp_reg = regionprops(lm, this.regprops); %gives the labels for the areas
                if ~isempty(temp_reg)
                    if (1 == 0) %debug plotting
                        bin_im = zeros(this.height, this.width);
                        %bin_im(temp_area.PixelIdxList) = 255;
                        figure; imshow(lm, [0 max(lm(:))]); colormap(jet);
                    end
                    [~, areai] = sort([temp_reg.Area], 'descend');
                    temp_reg = temp_reg(areai); %reorder in terms of area
                    area_sel = ([temp_reg.Area] >= this.MIN_SIZE_THRESH & [temp_reg.Area] <= this.MAX_SIZE_THRESH);
                    temp_reg= temp_reg(area_sel);
                    for kk=1:length(temp_reg)
                        temp_reg(kk).Orientation = -temp_reg(kk).Orientation; %the returned orientation is not correct, need to correct
                    end
                    this.nblobs(fi) = min(this.maxBlobs, length(temp_reg)); %take the biggest N blobs
                    this.areas(segFrames(ii),1:this.nblobs(fi)) = temp_reg(1:this.nblobs(fi)); %only keep the top sized blobs
                    this.assignBlobIDs(segFrames(ii)); % give the blobs unique IDs
                    this.detectTail(fi);
                    this.bodyCOM(fi,:) = this.computeBodyPos(fi,1); % 1) get the center of mass of animal
                    [this.nosePos(fi,:), this.noseblob(fi)] = this.findNose(fi); % 2) get the nose blob
                   
                    if fi>1 %compute frame-by-frame velocities
                        this.bodyVel(fi, :) = this.bodyCOM(fi,:) - this.bodyCOM(fi,:);
                        this.noseVel(fi, :) = sqrt(sum((this.nosePos(fi,:)-this.nosePos(fi-1,:)).^2));
                        if isnan(this.nosePos(fi-1,1))
                            vel = [0 0];
                        else
                            vel = this.noseVel(fi, :);
                        end
                    else
                        vel = [0 0];
                    end
                    
                    if dbg %debug plotting
                        bin_im = zeros(this.height, this.width);
                        for kk=1:length(temp_reg)
                            bin_im(temp_reg(kk).PixelIdxList) = kk;
                        end
                        figure; imshow(label2rgb(bin_im,'jet','k'));
                        hold on; plot(this.COM(fi,1,this.tailblob(fi)), this.COM(fi,2,this.tailblob(fi)),'r+');
                        hold on; plot(this.bodyCOM(fi,1), this.bodyCOM(fi,2),'yx','MarkerSize',12);
                    end
                else %this is a kluge - unclear what the best thing to do is if you don't detect a blob
%                     subI = max(1, ii-1);
%                     this.areas(segFrames(ii)) = this.areas(segFrames(subI)); %this will error on ii=1
%                     if fi < this.nFrames
%                         this.kf.s(fi+1) = kalmanf(this.kf.s(fi)); %let's update the Kalman Filter anyway
%                     end
                end
            end
        end
        
        
        
        % ------------------------------------------------------------------------------------------------------
        function newPos = correctDetection(this, frame, detectedPos, predictedPos)
            err_thresh = 7;
            areas = this.areas(frame,:);
            centers = [areas.Centroid];
            centers = reshape(centers, 2,[])';
            D = pdist2(predictedPos, centers); %distances from prediction to all centers
            [mind, mini] = nanmin(D);
            prevPos = [0 0];
            if frame > 1
                prevPos = this.nosePos(frame-1,:);
            end
            if isnan(this.noseblob(frame)) && mind < err_thresh %if we missed the nose blob for some reason
                newPos = centers(mini, :);
                this.nosePos(frame, :) = newPos;
                this.noseblob(frame) = mini;
                this.kf.s(frame).z = [newPos newPos-prevPos]';
                disp(sprintf('Corrected frame %i from %d, %d   to %d, %d', frame, detectedPos(1), detectedPos(2), newPos(1), newPos(2)));
            elseif mini ~= this.noseblob(frame) && mind < err_thresh %if a different blob is closer to the prediction
                % then we will call that the nose instead
                newPos = centers(mini, :);
                this.nosePos(frame, :) = newPos;
                this.noseblob(frame) = mini;
                this.kf.s(frame).z = [newPos newPos-prevPos]';
                disp(sprintf('Corrected frame %i from %d, %d   to %d, %d', frame, detectedPos(1), detectedPos(2), newPos(1), newPos(2)));
            elseif mini == this.noseblob(frame) && mind > err_thresh
                % do nothing for now
                % eventually want to rethreshold the frame
                
            end
            
        end

        % ------------------------------------------------------------------------------------------------------
        function propogateNosePosition(this, seedFrame)
        % function propogateNosePosition(this, seedFrame)
        %
        % The assigns the nose position to be the area with the same blob ID above and below the seed
        % frame
            noseID = this.blobID(seedFrame, this.noseblob(seedFrame));
            currFrame = seedFrame + 1;
            exit = 0;
            while (currFrame <= this.nFrames) && ~exit %forwards
                ids = this.blobID(currFrame,:);
                matchi = find(ids == noseID);
                if ~isempty(matchi)
                    this.noseblob(currFrame) = matchi;
                    this.nosePos(currFrame,:) = this.areas(currFrame, matchi).Centroid;
                    currFrame = currFrame + 1;
                else
                    exit = 1;
                end
            end
            currFrame = seedFrame -1;
            exit = 0;
            while (currFrame > 0) && ~exit %backwards
                ids = this.blobID(currFrame,:);
                matchi = find(ids == noseID);
                if ~isempty(matchi)
                    this.noseblob(currFrame) = matchi;
                    this.nosePos(currFrame,:) = this.areas(currFrame, matchi).Centroid;
                    currFrame = currFrame - 1;
                else
                    exit = 1;
                end
            end  
        end
        
        % ------------------------------------------------------------------------------------------------------
        function refineTracking(this, frames)
        % function refineTracking(this, frames)
        %
        % It makes sense to group these because they should all be performed together, and recomputed if the obh.
            this.removeStaticObjects();
            this.detectTail(frames);
            this.bodyCOM(frames,:) = this.computeBodyPos(frames,0);
            this.bodyOrient(frames) = this.computeBodyOrientation(frames);
            this.bodyOrient(frames) = this.fixOrientation(frames); %this should straighten out the orientations
            this.nosePos(frames, :) = this.findNose(frames); 
            this.vel(frames, :) = this.computeVelocity(frames);
        end
        
        % -----------------------------------------------------------------------------------------------------
        function detectTail(this, frames)
            % Trying to detect if we can see a tail on a frame by frame basis
            for ii=1:length(frames)
                fi = frames(ii);
                regions = this.areas(fi,:);
                maxo = 0;
                if (fi>1) && this.tailVisible(fi-1) %finding similarity with previously id'd tail 
                    prev_tail = this.areas(fi-1,this.tailblob(fi-1));
                    overlap = regionOverlap(prev_tail, regions, this.regprops);
                    [maxo, maxoi] = max(overlap);
                end
                if maxo > .1 %if similar, we've found another tail - this is a conservative threshold. Comparing the 
                    % overlap of id'd tail blobs to non-tail blobs, they are all down below
                    % .05...essentially zero.  This is really conservative - assuming that you get the
                    % tail in the first place.
                    this.tailVisible(fi) = 1;
                    this.tailblob(fi) = maxoi;
                else
                    % going to find the blob with the highest eccentricity to be the tail
                    ecc = [regions.MajorAxisLength]./[regions.MinorAxisLength];
                    [max_ecc, maxi] = nanmax(ecc);
                    if max_ecc>4.5 && regions(maxi).Area > 60%this is just an empirically defined threshold value
                        this.tailVisible(fi) = 1;
                        this.tailblob(fi) = maxi;
                    else
                        this.tailVisible(fi) = 0;
                        this.tailblob(fi) = NaN;
                    end
                end
            end
        end
        % ------------------------------------------------------------------------------------------------------
        function [nosePos, noseblob] = findNose(this, frames)
        % function nosePos = findNose(this, frames, blobi)
        %
        % So, the method for finding which blob is the nose is to take the tail and body center and take the 
        % distance along that vector for every blob centroid. The one with the largest is the nose.
            %nosePos = NaN* zeros(length(frames), 2);
            bestFrames = false(size(frames));
            candidate = zeros(size(frames));
            noseDists = zeros(size(frames));
            nosePos = this.nosePos(frames,:);
            noseblob = this.noseblob(frames);
            dot_thresh = 2000; % This is just an empirical value for the dot product of body vect/nose vector
            % go through the frame set and identify the best candidates for the nose
            for ii = 1:length(frames)
                fi = frames(ii);
                if this.tailVisible(fi) && (this.nblobs(fi) > 2)%we can only do this first method if the tail is present
                    % compute the vector from tail to body center
                    bodyVect = this.bodyCOM(fi,:) - this.areas(fi,this.tailblob(fi)).Centroid;
                    dist = [];
                    for jj=1:this.nblobs(fi)
                        areaVect = this.areas(fi,jj).Centroid - this.areas(fi,this.tailblob(fi)).Centroid;
                        dist(jj) = dot(areaVect, bodyVect);
                    end
                    [maxd, maxi] = max(dist); %the maximum distance and ind along tail-body vector of each blob
                    candidate(ii) = maxi;
                    noseDists(ii) = maxd;
                    bestFrames(ii) = maxd >= dot_thresh;
                end
            end
            bestFrames = find(bestFrames == 1);
            %bestFrames = frames(bestFramesLocal); % necessary to get the right overall frame numbers if called with a subset of frames
            for ii = 1:length(bestFrames)
                fi = bestFrames(ii);
                nosePos(fi,:) = this.areas(frames(fi), candidate(fi)).Centroid;
                noseblob(fi) = candidate(fi);
                this.nosePos(frames(fi),:) = nosePos(fi,:);
                this.noseblob(frames(fi)) = noseblob(fi);
                this.propogateNosePosition(frames(fi));
            end
 
            %let's assign the vectors
            this.nosePos(frames,:) = nosePos;
            this.noseblob(frames) = noseblob;
        end
        
        % ------------------------------------------------------------------------------------------------------
        function bodyCOM = computeBodyPos(this, frames, includeTail)
            %
            if isempty(includeTail) includeTail=0; end
            bodyCOM = NaN*zeros(length(frames), 2);
            for i = 1:length(frames)
                fi = frames(i);
                taili = this.tailblob(fi);
                if (isnan(taili)) taili = this.maxBlobs + 1; end %just assign it out of range for check
                temp_pos = [];
                temp_areas = this.areas(fi,1:this.nblobs(fi));
                for j=1:this.nblobs(fi)
                    this.orient(fi, j) = [this.areas(fi,j).Orientation]./ 180 * pi; %this is the rough estimate
                    positions = permute(reshape([temp_areas(j).Centroid], 2,[]), [3 1 2]);
                    this.COM(fi, :, j) = positions;
                    if (includeTail) || (j ~= taili) %exclude the tail in the body position calculation
                        %temp_pos = cat(1, temp_pos, temp_areas(j).PixelList);
                        temp_pos = cat(1, temp_pos, positions);
                        
                        %this.orient(fi, j) = [this.areas(fi,j).Orientation]./ 180 * pi; %this is the rough estimate
                    end
                end
                %bodyCOM(i, :) = mean(this.COM(fi,:,:),3);
                if isempty(temp_pos)
                    bodyCOM(i,:) = [NaN NaN];
                else
                    bodyCOM(i, :) = mean(temp_pos,1);
                end
            end
        end
        
       
        % ------------------------------------------------------------------------------------------------------
        function orient = computeBodyOrientation(this, frames)
        % function [this, orient] = computeBodyOrientation(this, frames)
        %
        % The idea here is to fill in the body space using dilation of the individual body blobs (excluding tail)
        % in order to get a reading of the animal orientation.
            se = strel('disk',20);
            orient = NaN*zeros(length(frames),1);
            for ii=1:length(frames)
                fi = frames(ii);
                if (this.nblobs(fi) > 1)
                    areas = this.areas(fi,1:this.nblobs(fi));
                    if (this.tailVisible(fi))
                        sel = true(this.nblobs(fi), 1);
                        sel(this.tailblob(fi)) = 0;
                        areas = areas(sel);
                    end
                    bw = false(this.height, this.width);
                    for jj=1:length(areas)
                        bw(areas(jj).PixelIdxList) = 1;
                    end
                    bw2 = imdilate(bw,se);
                    propstr = {'Orientation'};
                    props = regionprops(bw2, propstr);
                    if (~isempty(props))
                        orient(ii) = props(1).Orientation ./ 180 .* pi;
                    else 
                        orient(ii) = NaN;
                    end
                end
            end
            %this.orient(frames) = orient;
        end
        
        % ------------------------------------------------------------------------------------------------------
        function vel = computeVelocity(this, frames)
            % function vel = getVelocity(this, frames)
            % Computes the velocity as the frame by frame difference in position
            
            %in order to get a velocity for the first frame (0,0) position is assumed at frame 0 
            %this.vel = NaN*zeros(size(this.bodyCOM,1),size(this.bodyCOM,3));
            if isempty(frames) frames=1:this.nFrames; end
            if frames(1) > 1 %compute velocity from the frame before the range to get all frames requested
                fr = [frames(1)-1; frames(:)];
            else
                fr = [1; frames(:)];
            end
            %pos = this.bodyCOM(fr,:);
            pos = this.nosePos(fr,:);
            diff_pos = diff(pos, 1,1);%differences in each x,y position
            vel = sqrt(sum(diff_pos.^2, 2));
            this.vel(fr(2:end),1) = vel; % vel is 1 shorter than the fr
            this.direction = cart2pol(diff_pos(:,1), diff_pos(:,2)); %this is the direction of motion
            this.noseVel = this.vel;
            
            pos = this.bodyCOM(fr,:);
            diff_pos = diff(pos, 1,1);
            vel = sqrt(sum(diff_pos.^2, 2));
            this.bodyVel(fr(2:end),1) = vel; % vel is 1 shorter than the fr
            %diff_com = [COM(1,:); diff_com]; % add the first position to get an equal sized vector
            %this.direction(frames) = cart2pol(diff_com(:,1), diff_com(:,2)); 
        end
        
        % ------------------------------------------------------------------------------------------------------
        function orientation = fixOrientation(this, frames)
            % What we need to do from the ellipse orientation is to figure out the head direction
            %diffi = find(abs(this.orient(frames,:) - this.direction(frames, :)) > pi/2);
            
            % The general approach is to use the presence of a tail or motion vector to disambiguate the 
            % orientation of the animal
            
            %for frames with visible tail
            tp = [this.tailVisible];
            tp = logical(tp);
            %tp(:) = 0; %let's avoid some issues for now
            %tp = (tp ==1); %make it into a logical array
            tp = tp(frames);
            if(sum(tp))
                ntf = length(frames(tp));
                ftp = frames(tp);
                tailCOM = zeros(ntf,2);
                for ii = 1:ntf
                    taili = this.tailblob(ftp(ii));
                    tailCOM(ii,:) = this.areas(ftp(ii), taili).Centroid;
                end
                tailV = tailCOM - this.bodyCOM(ftp,:); %vector towards the tail
                [theta, ~] = cart2pol(tailV(:,1), tailV(:,2));
                %we want these to be close, the opposite of the tail vector and the orientation
                od = circularDistance(this.bodyOrient(ftp), theta+pi); 
                dif = abs(od) > pi/2;  % these estimates are way off, so rotate them 180deg
                this.bodyOrient(ftp(dif)) = mod(this.bodyOrient(ftp(dif)) + pi, 2*pi);
                switchedframes = ftp(dif);
                for ii=1:length(switchedframes) %propogate that change
                    startf = switchedframes(ii);
                    frame = min(startf+1, this.nFrames);
                    odf = abs(circularDistance(this.bodyOrient(startf), this.bodyOrient(frame)));
                    while(odf > pi/2)
                        this.bodyOrient(frame) = mod(this.bodyOrient(frame) + pi, 2*pi); %rotate 180deg
                        prev = frame;
                        frame = min(frame+1, this.nFrames);
                        odf = abs(circularDistance(this.bodyOrient(prev), this.bodyOrient(frame)));
                    end
                    % now propogate change in the opposite direction
                    startf = switchedframes(ii);
                    frame = max(1, startf-1);
                    odf = abs(circularDistance(this.bodyOrient(startf), this.bodyOrient(frame)));
                    while(odf > pi/2)
                        this.bodyOrient(frame) = mod(this.bodyOrient(frame) + pi, 2*pi); %rotate 180deg
                        prev = frame;
                        frame = max(frame-1, 1);
                          odf = abs(circularDistance(this.bodyOrient(prev), this.bodyOrient(frame)));
                    end
                end

                vel_thresh = nanmax(this.vel(2:end)/2);
                vel_thresh = 200;
                tp = this.vel(2:end) > vel_thresh;  % %threshold in px/sec
                tp = logical([0; tp(:)]);
                tp = tp(frames);
                ntf = length(frames(tp));
                ftp = frames(tp);
                vel_v = this.direction(ftp);
                od = circularDistance(this.bodyOrient(ftp), vel_v); %we want these to be close
                dif = abs(od) > pi/2;  % these estimates are way off, so rotate them 180deg
                this.bodyOrient(ftp(dif)) = mod(this.bodyOrient(ftp(dif)) + pi, 2*pi);
                switchedframes = ftp(dif);
                for ii=1:length(switchedframes) %propogate that change
                    startf = switchedframes(ii);
                    frame = min(startf+1, this.nFrames);
                    odf = abs(circularDistance(this.bodyOrient(startf), this.bodyOrient(frame)));
                    while(odf > pi/2)
                        this.bodyOrient(frame) = mod(this.bodyOrient(frame) + pi, 2*pi); %rotate 180deg
                        prev = frame;
                        frame = min(frame+1, this.nFrames);
                        odf = abs(circularDistance(this.bodyOrient(prev), this.bodyOrient(frame)));
                    end
                    % now propogate change in the opposite direction
                    startf = switchedframes(ii);
                    frame = max(1, startf-1);
                    odf = abs(circularDistance(this.bodyOrient(startf), this.bodyOrient(frame)));
                    while(odf > pi/2)
                        this.bodyOrient(frame) = mod(this.bodyOrient(frame) + pi, 2*pi); %rotate 180deg
                        prev = frame;
                        frame = max(frame-1, 1);
                          odf = abs(circularDistance(this.bodyOrient(prev), this.bodyOrient(frame)));
                    end
                end
            end
            orientation = this.bodyOrient(frames);
        end
        
        % ------------------------------------------------------------------------------------------------------
        function frames = timesToFrames(this, time_range)
            % This is a utility function that returns a column vector of frame numbers for a set of time ranges,
            % specified by an n by 2 matrix of times
            if isempty(time_range) %all times
                time_range = [this.times(1) this.times(end)];
            end
            frames = [];
            for ii = 1:size(time_range,1) %the number of rows in time_range
                fi = find(this.times >= time_range(ii,1) & this.times <= time_range(ii,2));
                frames = cat(1, frames, fi(:));
            end
        end
        
        % ------------------------------------------------------------------------------------------------------
        function trimFields(this)
            % this function shortens the properties fields to nFrames
            this.COM = this.COM(1:this.nFrames,:,:);
            this.bodyCOM = this.bodyCOM(1:this.nFrames,:,:);
            this.noseblob = this.noseblob(1:this.nFrames,:);
            this.tailblob = this.tailblob(1:this.nFrames,:);
            this.tailVisible = this.tailVisible(1:this.nFrames);
            this.bodyOrient= this.bodyOrient(1:this.nFrames,:);
            this.orient= this.orient(1:this.nFrames,:,:);
            this.vel = this.vel(1:this.nFrames,:);
            this.direction= this.direction(1:this.nFrames,:);
            this.times = this.times(1:this.nFrames);
            this.areas = this.areas(1:this.nFrames,:);
            this.nblobs = this.nblobs(1:this.nFrames);
            
        end
        
        % ------------------------------------------------------------------------------------------------------
        function trimMovie(this, frames)
            % this function shortens the tracked segment and all object data to only those
            % frames specified. DO NOT USE TO SUBSAMPLE FRAMES.  WILL BREAK THINGS
            frames = intersect(1:this.nFrames, frames); %just so that some bad call doesn't f things up
            this.nFrames = length(frames);
            this.frameRange = this.frameRange(frames);
            
            this.COM = this.COM(frames,:,:);
            this.bodyCOM = this.bodyCOM(frames,:);
            this.noseblob = this.noseblob(frames,:);
            this.nosePos = this.nosePos(frames,:);
            this.tailblob = this.tailblob(frames,:);
            this.tailVisible = this.tailVisible(frames);
            this.orient= this.orient(frames,:,:);
            this.bodyOrient = this.bodyOrient(frames,:);
            this.vel = this.vel(frames,:);
            this.direction= this.direction(frames,:);
            this.times = this.times(frames);
            this.areas = this.areas(frames,:);
            this.nblobs = this.nblobs(frames);
            
        end
        
        
        % ------------------------------------------------------------------------------------------------------
        function removeStaticObjects(this)
            % function removeStaticObjects(this)
            % 
            % Goes through the areas field of detected objects and removes
            % those that are similarly present for a significant number of
            % frames.
            limit = round(this.frameRate * 3); %X sec; number of frames necessary to be eliminated as static
            dist_thresh = 1; %What do we call the same position, along each axis 
            area_thresh = 15; %the amount that an area can change frame to frame to be counted as the same
            abs_area_thresh = 30; %the upper bound for the size for deleted objects
            max_block_size = 2^11; %the maximum number of frames processed at once - need to break up for long segments
                                   %for memory reasons
            nblocks = ceil(this.nFrames/max_block_size);
            this.blobsToDelete = []; di = 0;
            for jj = 1:nblocks
                if jj == nblocks
                   sz = this.nFrames - (nblocks-1)*max_block_size;
                else
                   sz = max_block_size;
                end
                bi = (jj-1)*max_block_size+1;
                frames = (bi:(bi+sz-1))';
                framem = repmat(frames, 1, this.maxBlobs); %this is the frame number for each element
                % First of all, go through and find the common positions
                all_pos = [this.areas(frames,:).Centroid];
                all_pos = reshape(all_pos, 2, sz, this.maxBlobs); 
                %all_pos = permute(all_pos, [2 3 1]);
                %all_pos = permute(all_pos, [1 3 2]);
                all_pos = reshape(all_pos, 2, []);
                all_pos = all_pos';%makes it an blob*frame x 2 matrix, in a all frames blob1, then all frames blob2, etc order
                %all_pos = reshape(all_pos, 2, []);
                %all_pos = reshape(all_pos, [], 2); 
                nz = logical(all_pos(:,1)) | logical(all_pos(:,2)); %logical for indexing sample/blob
                nzi = find(nz);
                nz_pos = [all_pos(nz,1) all_pos(nz,2)];
                distm = ipdm(single(nz_pos));
                distm = distm + (dist_thresh+1)*eye(size(distm,1)); %makes the diagonal fall above thresh so that they aren't zero.
                dist_sel = (distm <= dist_thresh); 
                dist_rc = find(sum(dist_sel,2) > limit);
                % let's also find an area difference to make sure they are the same blob
                all_area = [this.areas(frames,:).Area];
                nz_area = all_area(nz);
                area_diff = bsxfun(@minus, nz_area, nz_area');
                area_diff = area_diff + eye(size(area_diff,1))*area_thresh; %again, offset the same area comparison
                area_sel = (abs(area_diff) < area_thresh);
                area_rc = find(sum(area_sel,2) > limit);
                %absolute area
                abs_area_rc = find(nz_area < abs_area_thresh);

                %Now, go through and find rows that are too populated
                sel = dist_sel & area_sel;
                rowcounts = sum(sel,2);

                static_rows = find(rowcounts > limit);
                static_rows2 = intersect(dist_rc, area_rc);
                static_rows2 = intersect(static_rows2, abs_area_rc);
                % now get the frames that we want to delete
                %origi = frames(nzi(static_rows2));
                for i=1:length(static_rows2)
                    origi = nzi(static_rows2(i));
                    framei = framem(origi);
                    blobi = floor((origi-1)/sz) + 1;
                    %framei = floor(origi/this.maxBlobs);
                    %blobi = origi - (framei*this.maxBlobs) + 1;
                    %blobi = floor((origi-1)/this.nFrames) + 1;
                    %framei = origi - ((blobi-1)*this.nFrames);
                    di = di+1;
                    this.blobsToDelete(di,:) = [framei, blobi];
                    %this.deleteArea(framei, blobi);
                end
            end
            for jj = 1:size(this.blobsToDelete,1)
                this.deleteArea(this.blobsToDelete(jj,1), this.blobsToDelete(jj,2));
            end
        end
        % =-------------------------------------------------------------------------------------------
        function deleteArea(this, framei, blobi)
           % function deleteArea(this, framei, blobi)
           % 
           % Just deletes the blob by saving an empty structre in its place
           props = struct('Area',0, 'Centroid',[0 0],'BoundingBox', [0 0 0 0], 'MajorAxisLength',0,...
                'MinorAxisLength',0, 'Orientation',0,'Extrema',[0 0 0 0], 'PixelIdxList',[], 'PixelList',[],'Perimeter', 0);
           this.areas(framei,blobi) = props; 
           this.nblobs(framei) = this.nblobs(framei)-1;
        end
        % ------------------------------------------------------------------------------------------------------
        function clearCalcData(this)
         % function clearCalcData(this)
         % 
         % This is a function that initializes or reintializes the computed data
            this.COM = NaN*zeros(this.nFrames, 2, this.maxBlobs);
            this.orient = NaN*zeros(this.nFrames, this.maxBlobs);
            this.vel = NaN*zeros(this.nFrames, 1);
            this.direction = NaN*zeros(this.nFrames, this.maxBlobs);
            this.tailVisible = zeros(this.nFrames, 1);
            this.nosePos = NaN*zeros(this.nFrames, 2); %nose position
            this.nblobs = zeros(this.nFrames, 1);
            this.tailblob = NaN*zeros(this.nFrames,1);
            this.noseblob = NaN*zeros(this.nFrames,1);
            this.bodyOrient = NaN*zeros(this.nFrames, 1);
            this.bodyCOM = NaN*zeros(this.nFrames, 2);
            this.exploredProp = zeros(this.nFrames, 2);
            %initialize areas - must get these in the correct order too
            tempareas = struct('Area',0, 'Centroid',[0 0],'BoundingBox', [0 0 0 0], 'MajorAxisLength',0,...
                'MinorAxisLength',0, 'Orientation',0,'Extrema',[0 0 0 0], 'PixelIdxList',[], 'PixelList',[],'Perimeter', 0);
            this.areas = repmat(tempareas, this.nFrames, this.maxBlobs); % make a struct arraay length of nFrames
            this.initKalmanFilter();
            this.blobID = NaN*zeros(this.nFrames, this.maxBlobs);
            this.blob_num = 1;
        end
        
        % ------------------------------------------------------------------------------------------------------
        function clearPathData(this)
         % function clearPathData(this)
         % Clears the path data
            this.paths = [];
        end
        % ------------------------------------------------------------------------------------------------------ 
        function detectPaths(this, time, absoluteTime, useAvgFrame, imAx)
            % Essentially just calls detectEdgesInFrame twice, once for
            % each path, and saves the results
            pb = 0;
            
            [this.paths, rew_image] = detectEdgesInFrame(this, time, absoluteTime, useAvgFrame, imAx);
            [this.paths(2), distract_image] = detectEdgesInFrame(this, time, absoluteTime, useAvgFrame, imAx);
            
            if (pb)
                order = [2 1 3];
                r = this.avgFrame; g = this.avgFrame; b = this.avgFrame;
                r(distract_image) = 255; g(distract_image) = 0; b(distract_image) = 0;
                r(rew_image) = 0; g(rew_image) = 255; b(rew_image) = 0;
                colorim = cat(3, r, g, b);
                figure;
                imshow(colorim);
                hold on;
            end
            
        end
         % ------------------------------------------------------------------------------------------------------ 
        function detectRefinePaths(this,time,absoluteTime,useAvgFrame, imAx)
            this.detectPaths(time, absoluteTime, useAvgFrame, imAx);
            this.refinePaths(1);
            this.refinePaths(2);
            this.generatePathVertices();
        end
        
         % ------------------------------------------------------------------------------------------------------ 
        function makePathsSkel(this)
            if isempty(this.fullPaths)
                this.fullPaths = this.paths;
            end
            for ii = 1:length(this.paths)
                this.paths(ii) = skeletonizePath(this.paths(ii), this.width, this.height);
            end
            this.generatePathVertices();
        end
         % ------------------------------------------------------------------------------------------------------ 
        function makePathsFull(this)
            if ~isempty(this.fullPaths)
                this.paths = this.fullPaths;
                this.fullPaths = [];
                this.generatePathVertices();
            else
                disp('There is no fullPaths field.');
            end
        end
        % ------------------------------------------------------------------------------------------------------
        function makeTrailShadow(this, shadowSize, pathNum)
            if ~isempty(this.paths)
                this.makePathsSkel();
                this.trailShadow = findAllPointsWithinDistance(this.paths(pathNum).PixelList, 20);
            end
        end
        % ------------------------------------------------------------------------------------------------------
        function generatePathVertices(this)
            % function generatePathAreas(this)
            %
            % So, detecting the distance of a floating point (effectively analog) position from a pixel
            % based area is not as simple as euclidean distance. Each pixel is defined by a "position" that 
            % the center of that square but is actually an area.  In order to overcome that, we use positional
            % vertices to define the area that we want to compute the distance from.  This computes that
            % area from a list of pixel positions.
            this.pathVertices = this.paths;
            for ii = 1:length(this.paths)
                px = this.paths(ii).PixelList;
                vert = NaN*zeros(size(px,1)*4, 2); %upper bound is 4 times the px
                vi = 1;
                for jj = 1:size(px,1)
                    pxt = px(jj,:);
                    vert(vi,:)   = [pxt(1)-.5, pxt(2)-.5]; %left, top
                    vert(vi+1,:) = [pxt(1)+.5, pxt(2)-.5]; %right, top
                    vert(vi+2,:) = [pxt(1)-.5, pxt(2)+.5]; %left, bottom
                    vert(vi+3,:) = [pxt(1)+.5, pxt(2)+.5]; %right, bottom
                    vi = vi+4;
                end
                uv = unique(vert, 'rows', 'R2012a');
                num = size(uv,1);
                uvx = sub2ind(size(this.avgFrame), uv(:,2), uv(:,1));
                this.pathVertices(ii).Area = num;
                this.pathVertices(ii).PixelList = uv;
                this.pathVertices(ii).PixelIdxList = uvx;
            end
            if ~isempty(this.pathVertices)
                this.pathVerticesDefined = 1;
            end
        end
        % ------------------------------------------------------------------------------------------------------
        function [e, eimage] = detectEdgesInFrame(this, time, absoluteTime, useAvgFrame, imAx)
            % function e = detectEdgesInFrame(this, time, absoluteTime)
            %
            % Detects the edges within an movie frame. This is the MouseTracker object, time is the time during the
            % movie and absoluteTime is a boolean flag to use movie relative time or absolute time (useful for detecting
            % edges in parts of the movie that aren't used for tracking).
            
            EDGE_LEN_THRESH = 20;
            disk_size = 20;
            pb = 0;
            % choose the image to use, then adjust to maximize contrast
            if (isempty(useAvgFrame))
                useAveFrame = 0;
            end
            if ~useAvgFrame
                if ~absoluteTime
                    f = this.timesToFrames([time time+1]);
                else % use the absolute time of the video
                    f = time*this.frameRate+1;
                    f = round(f)-this.frameRange(1);
                    f = [f f+1];
                end
                vid_struct = this.readFrames(f);
                gf = vid_struct.frames(1).cdata;
            else
                gf = this.avgFrame;
            end
            mingf = double(min(gf(:))); maxgf = double(max(gf(:)));
            gf = imadjust(gf, [mingf/255; maxgf/255], [0; 1]);
            
            % get user input about which lines they want to be marked
            if isempty(imAx)
                fh = figure; imshow(gf);
                [cx, cy] = ginput; %get 2 points from user interaction
                close(fh);
            else
                axes(imAx); imshow(gf);
                [cx, cy] = ginput; 
            end
            
            [ei, thresh] = edge(gf, 'canny'); %first detect edges
            %thresh(1) = .8*thresh(2); %inefficient, to detect twice, but I want to adjust the thresh
            %ei = edge(gf, 'canny', thresh); %first detect edges
            ei = imclose(ei, strel('square', 3)); %morphological close, fills in small (1 px) gaps
            props = {'Area', 'PixelIdxList', 'PixelList'};
            e = regionprops(ei, props); %gets the binary regions defined by the edges
            % eliminate the short edges
            ea = [e.Area];
            long = ea >= EDGE_LEN_THRESH; 
            e = e(long);
            ea = ea(long);
            [~, order] = sort(ea); e = e(order); %sort by area
            
            %find areas where we've clicked
            p = round([cx(:), cy(:)]); %clicked points
            match = zeros(length(e),1); %areas that have been matched
            for ii=1:size(p,1) 
                for jj = 1:length(e)
                    linepx = e(jj).PixelList;
                    dists = sqrt((linepx(:,1) - p(ii,1)).^2 + (linepx(:,2) - p(ii,2)).^2); %distances to every point in area
                    if min(dists) <= 10 %then we've clicked very near to some edge area
                        match(jj) =1;
                    end
                end
            end
            e = e(match == 1); % only keep areas that have been matched
            e = mergeAreas(e); % want this user interaction to result in a single area
            % This is the returned binary image
            eimage = false(size(ei));
            eimage(e.PixelIdxList) = 1;
%             eimage = imfill(eimage, 'holes'); %We want to fill in the holes in the detected path
%             e = regionprops(eimage,props); %now just redetect.
%             e = mergeAreas(e);
            
            if (pb)
                hold on;
                overlay = gf;
                neg = gf;
                overlay(e.PixelIdxList) = 255;
                neg(e.PixelIdxList) = 0;
                colorim = cat(3, overlay, neg, neg);
                imshow(colorim);
            end
        end
        % ------------------------------------------------------------------------------------------------------
        function pathIm = plotPaths(this)
            %returns a binary image (uint8 format) of the detected paths
            
            pathIm = ones(this.height, this.width, 'uint8');
            for ii = 1:length(this.paths)
                path = this.paths(ii);
                pathIm(path.PixelIdxList) = 0;
                %pathIm = imfill(pathIm,'holes');
            end
        end
        % ------------------------------------------------------------------------------------------------------
        function pathIm = plotPathsOnBg(this, pathNums)
            if (~exist('pathNums','var'))  
                pathNums = 1:length(this.paths);
            elseif(~isempty(pathNums)) %if we're given a set of paths, use that, otherwise plot all of them
                pathNums = intersect(1:length(this.paths), pathNums); %make sure we don't try to plot anything not there
            else
                pathNums = 1:length(this.paths);
            end
            pathIm = cat(3, this.avgFrame, this.avgFrame, this.avgFrame);
            %pathIm = zeros(this.height, this.width, 3);
            color_order = [2 1 3];
            for ii = 1:length(pathNums)
                %set color layer to 255
                pi = pathNums(ii); 
                c = color_order(mod(pi-1, 3)+1);
                layer = pathIm(:,:,c);
                path = this.paths(pi);
                %path = skeletonizePath(path, this.width, this.height);
                layer(path.PixelIdxList) = 200;
                pathIm(:,:,c) = layer;
                oc = find(1:3 ~= c);
                %set the other layervidss to 0
                layer = pathIm(:,:,oc(1));
                layer(path.PixelIdxList) = 0;
                pathIm(:,:,oc(1)) = layer;
                layer = pathIm(:,:,oc(2));
                layer(path.PixelIdxList) = 0;
                pathIm(:,:,oc(2)) = layer;
            end
        end
        % ------------------------------------------------------------------------------------------------------
        function refinePaths(this, pathNum)
            % function removePaths(this)
            %
            % This pops up the background with paths plotted on it (the path specified)
            % so that the user can remove portions of it.
            props = {'Area', 'PixelIdxList', 'PixelList'};
            fh = figure;
            set(fh, 'WindowKeyPressFcn', @this.exitMovieLoop);
            this.exitMovie = 0;
            while(~this.exitMovie)
                imshow(this.plotPathsOnBg(pathNum));
                del_poly = impoly(gca);
                del_roi = del_poly.createMask();
                del_i = find(del_roi);
                path_i = this.paths(pathNum).PixelIdxList;
                [keep, ki] = setdiff(path_i, del_i);
                this.paths(pathNum).PixelList = this.paths(pathNum).PixelList(ki,:);
                this.paths(pathNum).PixelIdxList = keep;
                this.paths(pathNum).Area = length(keep);   
            end 
            eimage = false(this.height, this.width);
            eimage(this.paths(pathNum).PixelIdxList) = 1;
            eimage = imclose(eimage,strel('disk',3)); 
            eimage = imfill(eimage, 'holes'); %We want to fill in the holes in the detected path
            e = regionprops(eimage,props); %now just redetect.
            this.paths(pathNum) = mergeAreas(e);
            imshow(this.plotPathsOnBg([]));
        end
        % ------------------------------------------------------------------------------------------------------
        function [noseDist, vects, ret_frames] = noseDistToTrail(this, frames, pathNum, varargin)
        % function noseDist = noseDistToTrail(this, frames, pathNum, varargin)
        % 
        % returns the distances for the nose position to the closest point of trail
        % for each of the frames specified.
        % frames - the frames to be considered
        % pathNum - path number, generally 1 for rewarded, 2 for unrewarded
        % varargin - 1) the threshold distance for which to return the frames where the distance is within.
        % Always returns a POSITIVE distance
            if isempty(frames)
                frames = 1:this.nFrames;
            end
            if ~this.pathVerticesDefined
                this.generatePathVertices();
            end
            if (length(this.paths) >= pathNum)
                nosePos = this.nosePos(frames,:);
                trailPos = this.pathVertices(pathNum).PixelList;
                distm = ipdm(single(nosePos), single(trailPos));
                %distm = ipdm(nosePos, trailPos, 'Subset', 'SmallestFew', 'limit', 10);
                % Must figure out if the animals' nose is over the trail or
                % not. Do this by getting the 4
                % closest trail points to each nose position and see if they encircle it.
                noseDist = NaN*zeros(size(nosePos,1),6); 
                closestTrailP = ones(size(nosePos,1),2,6);
                for ii = 1:6
                    [noseDist(:,ii), mini] = nanmin(distm, [], 2); %get the minimum value
                    li = sub2ind(size(distm), (1:size(nosePos,1))', mini);
                    distm(li) = NaN; %set the mins to NaN to get the next closest on next iteration
                    closestTrailP(:,:,ii) = trailPos(mini, :);
                end
                over = isContained(nosePos, closestTrailP);
                
                % When the nose is close, the angles to whole number
                % vertices deviate from the continuous position of the nose
                % and the trail. Need to correct for this. When the
                % distance is <6px then the error is <5deg, which is
                % acceptable and under the normal turning variation.
%                 closei = find(noseDist < 6);
%                 for ii = 1:length(closei)
%                     up = ceil(nosePos);
%                     down = floor(nosePos);
%                     dists
%                     for jj = 1:2
%                     end
%                 end
                noseDist = noseDist(:,1);
                vects = closestTrailP(:,:,1) - nosePos;
                
                
                noseDist(over) = 0; %set those points to zero because the nose IS over the trail itself
            else % if there are no paths return zeros, but if there are return a mock path result
                if isempty(this.pathVertices)
                    noseDist = NaN*zeros(length(frames), 1);
                    vects = NaN*zeros(length(frames), 2);
                    disp([this.videoFN ': noseDistToTrail']);
                    disp('There must be a detected path to calculate distances');
                else %make a fake trail
                    nosePos = this.nosePos(frames,:);
                    fakeTrail = this.makeMirrorPath(this.pathVertices(1));
                    trailPos = fakeTrail.PixelList;
                    distm = ipdm(single(nosePos), single(trailPos));
                    [noseDist, mini] = nanmin(distm, [], 2);
                    closestTrailP = trailPos(mini, :);
                    vects = closestTrailP - nosePos;
                end
            end
            if ~isempty(varargin)
                thresh_dist = varargin{1}; %varargin 1 is the threshold distance to only return those distances under that
                %noseDist = noseDist * this.mm_conv;
                seli = abs(noseDist) < thresh_dist;
                noseDist = noseDist(seli) ;
                ret_frames = frames(seli);
                vects = vects(seli,:);
            else
                ret_frames = frames;
            end
        end
        % ------------------------------------------------------------------------------------------------------
        function vects = noseVectorToTrail(this, frames, pathNum)
        %function vects = noseVectorToTrail(this, frames, pathNum)   
        %
        % returns the vectors for each frame. This function computes this
        % vector separately for each frame, so it may be slow if used for
        % large numbers of frames at once.  Thus it is recommended that the
        % frames be selected prior to calling this method rather than
        % performing this on a whole movie worth of frames.
            if ~this.pathVerticesDefined
                this.generatePathVertices();
            end
            if (length(this.paths) >= pathNum)
                nosePos = this.nosePos(frames,:);
                trailPos = this.pathVertices(pathNum).PixelList;
                for ii = 1:size(nosePos,1)
                    distm = ipdm(nosePos, trailPos);
                end
            else
                vects = NaN*zeros(length(frames), 2);
                disp([this.videoFN ' : noseVectorToTrail']);
                disp('There must be a detected path to calculate distances');
            end
            
        
        end
        % ---------------------------------------------------------------------------------------------------
        function headingVect = headingFromBodyToNose(this, frames)
        %function headingVect = headingFromBodyToNose(this, frames)   
        %
            headingVect = this.nosePos(frames,:) -  this.bodyCOM(frames,:);
            theta = cart2pol(headingVect(:,1), headingVect(:,2));
            %this.bodyOrient(frames) = theta;
            headingVect = theta;
        end
        
        % ---------------------------------------------------------------------------------------------------
        function headingVect = headingFromMotion(this, frames)
        %function headingVect = headingFromBodyToNose(this, frames)   
        %
            if isempty(frames)
                frames = 1:this.nFrames;
            end
            headingVect = gaussianFilter(this.direction, 1,'conv');
            headingVect = headingVect(frames);
        end
        
        % ---------------------------------------------------------------------------------------------------
        function [fromOrthogonal, fromParallel, orthoTheta] = headingFromMotion_TrailRelative(this, frames, pathNum)
        %function [fromOthogonal, fromParallel] = headingFromMotion_TrailRelative(this, frames, pathNum)
        %
        % Two functional outputs:
        % 1) fromOthogonal - the angle relative to the orthogonal vector to
        % the closest point on the trail.
        % 2) fromParallel - the angle relative to the vector parallel to
        % the trail (orthogonal to the orthogonal vector). Importantly,
        % this vector is sign flipped based on the side of the trail the
        % animal is on, so that when the animal crosses the trail it stays
        % in the same direction rather than turning 180 deg.  
            pb=0;
            if isempty(frames)
                frames = 1:this.nFrames;
            end
            absHeading = this.headingFromMotion(frames);
            [noseDist, ~, orthoTheta] = this.orthogonalDistFromTrail(frames,pathNum);
            %orthoTheta = mod(cart2pol(vects(:,1), vects(:,2)),2*pi); %the orientation of the vector to the closest trail point
            fromOrthogonal = circ_dist(absHeading, orthoTheta);
            %[~, fromOrthogonal2] = rotationDirection(orthoTheta, absHeading);
            absCheck = mod(fromOrthogonal + orthoTheta, 2*pi);
            isright = (noseDist > 0); paraVect = orthoTheta;
            paraVect(isright) = orthoTheta(isright) + pi/2; paraVect(~isright) = orthoTheta(~isright) + 3*pi/2;
            fromParallel = circ_dist(absHeading, paraVect); 
            %headingVect = relToTrail;
            if pb
                this.plotNosePosition(frames);
                % Vector to subtract
                %[upara, vpara] = pol2cart(paraVect,abs(noseDist));
                %quiver(this.nosePos(frames,1), this.nosePos(frames,2), upara, vpara,'b', 'AutoScale', 'off');
                %hold on;
                [u, v] = pol2cart(orthoTheta,abs(noseDist));
                quiver(this.nosePos(frames,1), this.nosePos(frames,2), u, v,'b', 'AutoScale', 'off');
                hold on;
                % Red - absolute heading
                [u,v] = pol2cart(absHeading, 1);
                quiver(this.nosePos(frames,1), this.nosePos(frames,2), u, v, 'r','AutoScale', 'off');
                [u, v] = pol2cart(absCheck,1);
                quiver(this.nosePos(frames,1), this.nosePos(frames,2), u, v, 'k', 'AutoScale', 'off');
                % Magenta - relative heading
                [u,v] = pol2cart(fromParallel, 2);
                quiver(this.nosePos(frames,1), this.nosePos(frames,2), u, v, 'm', 'AutoScale', 'off');
                %[u,v] = pol2cart(fromOrthogonal2, 3);
                %quiver(this.nosePos(frames,1), this.nosePos(frames,2), u, v, 'c', 'AutoScale', 'off');
                scalea = linspace(0,2*pi, 9);
                scaler = (1:8) * 3 + 20;
                [u, v]  = pol2cart(scalea(1:end-1), scaler);
                quiver(50*ones(1,8), 50*ones(1,8), u, v, 'AutoScale', 'off');
            end
        end
        
        % ------------------------------------------------------------------------------------------------------
        function [dists, thetaDiff, orthoTheta] = orthogonalDistFromTrail(this, frames, pathNum)
        % function dists = orthogonalDistFromTrail(this, frames)
        % 
        % This function returns a signed value - leftward is negative, 
        % rightward is positive! Let's put several pieces together
            if isempty(frames) frames = 1:this.nFrames; end
            [noseDist, vects, ~] = this.noseDistToTrail(frames, pathNum);
            orthoTheta = mod(cart2pol(vects(:,1), vects(:,2)),2*pi); %the orientation of the vector to the closest trail point
            headingTheta = mod(this.headingFromBodyToNose(frames),2*pi); %the orientation of the animal
            [rotations, thetaDiff] = rotationDirection(headingTheta, orthoTheta); 
            dists = noseDist;
            dists = dists .* rotations; %gives it a sign
            
        end
        
        % ------------------------------------------------------------------------------------------------------
        function [crossings, dir, dists] = findTrailCrossings(this, frames, pathNum, varargin)
        % function [crossings, dir] = findTrailCrossings(this, frames, pathNum)
        %
        % This is the best function in the world. Finds trail crossings, and tells you the direction of
        % each.  From left->right is positive, from right->left is negative. 
        % Varargin provides additional functionality if you want to oversample the signal to get a more
        % precise time of crossing. t is that time.
        %
        % THIS MAY NOT WORK AS EXPECTED IF NOT CALLED WITH A CONTINUOUS SET OF FRAMES
            dist_thresh = 5;
            dists = this.orthogonalDistFromTrail(frames, pathNum);
            s1 = dists(1:end-1);
            s2 = dists(2:end);
            c = s1.*s2; % find the points where the distance from the trail changes sign, smart
            crossings = find(c < 0); 
            crossing_dists = dists(crossings); %also have to select the actual crossings, rather than the
                                               %mouse turning far from the trail.
            ci = abs(crossing_dists) < dist_thresh;
            crossings = crossings(ci);
            dir = dists(crossings) < 0;
            dists = dists(crossings);
            
            % remove the duplicates that may arise from very close on the
            % trail behavior (only count first crossing if there are
            % multiple crossings close by and near the trail
            if ~isempty(crossings)
                sep = diff(crossings);
                sep = [100; sep];
                sel = sep > 1;
                crossings = crossings(sel);
                dir = dir(sel);
                dists = dists(sel);
            end
        end    
        
        
        % -----------------------------------------------------------------------------------------------------
        function [turnFrames, dir, dists] = findFollowingTurns(this, frames, pathNum, threshDist, wind, varargin)
        % function [turnFrames, dir, dists] = findFollowingTurns(this, frames, pathNum, threshDist, wind)
        %
        % turnFrames - the list of frames where turns were detected
        % dir - direction of the turn. -1: rightward, 1: leftward
        % dists - the orthogonal distance from the trail in a window around
        % the turn.

            if nargin > 5
                pb = varargin{1};
            else
                pb = 0;
            end
            vel_conv = this.mm_conv * this.frameRate;
            if pb
                figure; hold on;
            end
            thresh_vel = [40 200];
            thresh_theta = [.1745 2.9671]; 
            filt_vel = gaussianFilter(this.noseVel, 2.5, 'conv') .* vel_conv;
            followingFrames = this.getFollowingSegments(frames, pathNum, threshDist);
            turnFrames = []; dir = []; dists = [];
            for ii = 1:size(followingFrames, 1)
                fr = followingFrames(ii,1):followingFrames(ii,2);
                vel_pass = filt_vel(fr) > thresh_vel(1) & filt_vel(fr) < thresh_vel(2);
                [dists, thetaDiff] = this.orthogonalDistFromTrail(fr, pathNum);
                theta_pass = abs(thetaDiff) >= thresh_theta(1) & abs(thetaDiff) <= thresh_theta(2);
                vel_pass = vel_pass & theta_pass;
                smooth_dists = gaussianFilter(dists, 2,'conv', 'valid'); %do the smoothing of position beforehand
                [~, turnInds,dirs] = findMinAndMaxPeaks(smooth_dists, -20, 2, 0, 0); %signal, threshold, window width, filter width
                vp_turns = vel_pass(turnInds);
                turnFrames = cat(1, turnFrames, turnInds(vp_turns) + fr(1) - 1);
                dir = cat(1, dir, dirs(vp_turns));
                %turns = turns & vel_pass; % just keep those turns that happen above a certain velocity
                %turnFrames = cat(1, turnFrames, find(turns) + fr(1) - 1);

                if pb        
                    subplot(2,1,1);
                    plot(this.times(fr), dists); 
                    plot(this.times(fr(turnInds)), dists(turnInds), 'ro');
                    hold on; plot(this.times(fr), smooth_dists, 'g'); 
                    xlabel('Time'); ylabel('Lateral distance from Trail'); 
                    
                    subplot(2,1,2);
                    hold on; plot(this.times(fr), this.headingFromMotion_TrailRelative(fr,pathNum));
                    xlabel('Time'); ylabel('Movement heading');
                end
            end 
            
            dists = NaN*zeros(length(wind), length(turnFrames)); %big preallocation
            for ii = 1:length(turnFrames)
                windi = wind + turnFrames(ii);
                inRange = find(windi >= 1 & windi <= this.nFrames);
                windi = windi(inRange);
                
                %centi = find(windi == turnFrames(ii), 1, 'first');
                %save_centi = find(wind == 0, 1, 'first');
                temp_dists = this.orthogonalDistFromTrail(windi, 1);
                dists(inRange, ii)=temp_dists;
            end
        end
        
        % -----------------------------------------------------------------------------------------------------
        function propNose = proportionFramesWithNose(this)
            np = this.nosePos(:,1);
            propNose = sum(~isnan(np))/this.nFrames;
        end
        
        % ------------------------------------------------------------------------------------------------------
        function meanDist = meanOrthogonalDistFromTrail(this, frames, pathNum)
        % function dists = meanOrthogonalDistFromTrail(this, frames, pathNum, threshDist)
        %c
        %
            dists = this.orthogonalDistFromTrail(frames, pathNum);
            meanDist = nanmean(dists);
        
        end
        
        % ------------------------------------------------------------------------------------------------------
        function dists = orthogonalDistFromTrailPerSection(this, frames, pathNum, threshDist)
        % function dists = meanOrthogonalDistFromTrail(this, frames, pathNum, threshDist)
        %
        %
            followingFrames = this.getFollowingSegments(frames, pathNum, threshDist);
            dists = zeros(size(followingFrames,1), 1);
            for ii = 1:size(followingFrames, 1)
                dists(ii) = this.meanOrthogonalDistFromTrail(followingFrames(ii,1):followingFrames(ii,2), pathNum);
            end
        end
        
        % ------------------------------------------------------------------------------------------------------
        function plotVectsToTrail(this, frames, pathNum)
            [noseDist, vects, ~] = this.noseDistToTrail(frames, pathNum);
            np = this.nosePos(frames,:);
            
            figure;
            imshow(this.plotPathsOnBg()); %make the background
            hold on;
            quiver(np(:,1), np(:,2), vects(:,1), vects(:,2), 0);
        end
        
            
        % ------------------------------------------------------------------------------------------------------
        function pTime = propTimeOnTrail(this, frames, trailNum, threshDist)
        %function pTime = propTimeOnTrail(this, frames, trailnum)    
        %
        % Return the percentage of the time period given that the mouse 
        % was within the threshold distance from the trail.
            %traillNum = 1; % Eventually we need 2+ trails
            if nargin < 4 %default distance
                threshDist = 10;
            end
            if isempty(frames)
                frames = 1:this.nFrames;
            end
            pTime = 0; %default return value
            dists = this.noseDistToTrail(frames, trailNum);
            nn = ~isnan(dists);
            nnframes = frames(nn);
            dists = dists(nn);
            numFrames = length(nnframes);
            closeFrames = sum(dists < threshDist);
            pTime = closeFrames/numFrames;
        end
        % ---------------------------------------------------------------------------------------------------
        function followingFrames = getFollowingSegments(this, frames, trailNum, threshDist)
        % function followingFrames = getFollowingSegments(this, frames, trailNum, threshDist)
        %
        % Returns the frame segments for which the mouse is following the trail
        
            if isempty(frames)
                frames = 1:this.nFrames;
            end
            %traillNum = 1; % Eventually we need 2+ trails
            if nargin < 4 %default distance
                threshDist = 20;
            end
            distTraveled = [];
            frameNums = [];
            
            dists = this.noseDistToTrail(frames, trailNum);
            nn = ~isnan(dists);
            nnframes = frames(nn);
            dists = dists(nn);
            numFrames = length(nnframes);
            trackingInds = find(dists < threshDist); %close to trail
            trackJumps = diff(trackingInds); %identify contiguous and non are 
            skips = find(trackJumps > 1); % noncontiguous
            if ~isempty(trackingInds)
                seq_end = [trackingInds(skips)' trackingInds(end)]'; %these are the indices of segments of following
                seq_start = [trackingInds(1) trackingInds(skips+1)']';
                frame_start = nnframes(seq_start); %now get the actual frame numbers
                frame_end = nnframes(seq_end);
                if (isempty(frame_end)) frame_end = nnframes(trackingInds(end)); end
                frameNums = [frame_start(:) frame_end(:)];
            end
            followingFrames = frameNums;
        end
        % ---------------------------------------------------------------------------------------------------
        function [distTraveled_ret, frameNums_ret] = distanceOnTrail(this, frames, trailNum, threshDist)
        % function [dists, frames] = distanceOnTrail(this, frames, threshDist) 
        %
        % Returns the distances travelled along the trail, along with the
        % frame numbers for the following onset and offset.
        distTraveled = [];
        frameNums = getFollowingSegments(this, frames, trailNum, threshDist);
        for ii=1:size(frameNums,1)
            range = frameNums(ii,1):frameNums(ii,2);
            pos = this.nosePos(range,:);
            pos = reshape(pos(~isnan(pos)),[],2); %this is to deal with the
            % fact that there may be nan position frames in the middle.
            point_diffs = diff(pos,1,1);
            if isempty(point_diffs)
                distTraveled(ii) = 0;
            else
                distTraveled(ii) = nansum(sqrt(nansum(point_diffs.^2,2)));
            end
        end
        if isempty(distTraveled)
            distTraveled = 0;
            frameNums = [1 1];
        end
        distTraveled_ret(ii) = distTraveled(ii);
        %distTraveled_ret = distTraveled * this.mm_conv;
        frameNums_ret = frameNums;
        end
        
        % ------------------------------------------------------------------------------------------------------
        function prop = propTrailFollowed(this, frames, trailNum, threshDist)
        % function propTrailFollowed(this, frames, trailNum, threshDist)
        %
        % This function returns the proporation of the trail pixels that the mouse came within a given
        % distance of.  Seems like a clean way of measuring the completeness of his trail exploration.
        % Algorthmically, this is a similar problem to finding the following segments except reversed.
        if isempty(frames)
            frames = 1:this.nFrames;
        end
        np = this.nosePos(frames,:);
        nn = ~isnan(np(:,1));
        np = np(nn,:);
        trailPos = this.paths(trailNum).PixelList; 
        distm = ipdm(single(np), single(trailPos));
        [trailDist, mini] = nanmin(distm, [], 1);
        explored = find(trailDist <= threshDist);
        npx = length(this.paths(trailNum).PixelIdxList);
        prop = length(explored)/npx;
        
        end
        % ------------------------------------------------------------------------------------------------------
        function turning_total = totalTurning(this, frames)
            measure_name = 'direction';
            measure = this.(measure_name);
            turning_vect = diff(measure);
            framei = frames > 1;
            frames = frames(framei);
            turning_vect = turning_vect(frames-1);
            turning_total = nansum(turning_vect);
        end
        
        % ------------------------------------------------------------------------------------------------------
        function nose_jumps = findNoseJumps(this, dist_thresh, frames)
            % reports the first frame after there is a large jump in the nose position
            includeNAN = 1;
            if isempty(frames)
                frames = 1:this.nFrames;
            end
            dists = sqrt(sum(diff(this.nosePos(frames,:)).^2, 2));
            nose_jumps = find(dists >= dist_thresh);
            if includeNAN
                ni = isnan(this.nosePos(frames,1));
            else
                ni = [];
            end
            nose_jumps = [nose_jumps; ni];
            nose_jumps = sort(unique(nose_jumps));
        end
        
        % ------------------------------------------------------------------------------------------------
        function showFrameAroundMouse(this, framei, movieType, varargin)
            wind_size = 150;
            np = this.nosePos(framei, :);
            
            if ~isnan(np(1))
                c = np;
            else
                cents = [];
                for ii = 1:nblobs(framei)
                    cents = cat(1, cents, this.areas(framei, ii).Centroid);
                end
                c = mean(cents);
            end
            dispCrop = [c(1) - wind_size/2,                  c(2) - wind_size/2, ...
                        this.width-( c(1)+wind_size/2 ),     this.height-( c(2)+wind_size/2 ) ];
            this.showFrame(framei, movieType, dispCrop, varargin{:});
        end
        
        % ------------------------------------------------------------------------------------------------------
        function showFrame(this, framei, movieType, dispCrop, varargin)
            % function showFrame(this, framei, useBinMovie)
            %
            % plots a frame of the movie,
            if isempty(dispCrop) dispCrop = [1 1 this.width this.height]; end
            if (nargin >= 5)
                logIm = varargin{1};
            else
                logIm = false(this.height, this.width);
            end
            dbg = 0;
            if strcmp(movieType, 'bin')
                bf = zeros(this.height, this.width, 'uint8');
                for jj=1:size(this.areas,2);
                    on = this.areas(framei,jj).PixelIdxList;
                    bf(on) = 1;
                end
                %also, need to highlight any areas to delete
                if dbg
                    fi = find(this.blobsToDelete(:,1) == framei);
                    for jj = 1:length(fi)
                        on = this.areas(this.blobsToDelete(fi(jj),1),this.blobsToDelete(fi(jj),2)).PixelIdxList;
                        bf(on) = bf(on)+jj;
                    end
                end
                % way to do it with white trails
                %pathIm = this.plotPaths*4;
                %bf = bf+pathIm + uint8(logIm)*3;
                %imshow(bf*255); hold on; 
                
                % way to do it with colored trails
                if isempty(this.colorBG)
                    this.colorBG = uint8(this.plotPathsOnBg()); %gives us a colored base
                end
                pathIm = this.colorBG;
                for ii = 1:3  
                    pathIm(:,:,ii) = (255*bf) + pathIm(:,:,ii);
                    if ii == 3 || ii==1
                        pathIm(:,:,ii) = uint8(logIm)*255 + pathIm(:,:,ii);
                    end
                end
                imshow(pathIm); hold on; 
                %imshow(label2rgb(bf, 'cool','k')); hold on; 
                
                %f = label2rgb(bf, 'cool','k');
                %pos = max(f,3); pos(logIm) = 255;
                %neg = max(f,3); f; neg(logIm) = 0;
            else
                f = this.readMovieSection(framei, movieType);
                f = this.plotPathsOnFrame(f, 1:2, 1);
                if (size(f,3) == 1)
                    f = this.gray2rgb(f);
                end
                if (nargin >=5) % A little overlay, hopefully
                    f(:,:,3) = 255 * uint8(logIm);
                    f(:,:,1) = 255 * uint8(logIm);
                end
                %neg = f; neg(logIm) = 0;
                %f = cat(3, neg, neg, pos);
                %f(:,:,3) = max(1, f(:,:,3)+cast(logIm,'uint8')*255); %not quite sure what I'm doing here
                %f(:,:,1) = 
                imshow(f);
                hold on;
            end
            %annotate the image
            if ~isempty(this.areas)
                hold on;
                for jj = 1:size(this.COM,3)
                    plot(this.COM(framei,1,jj), this.COM(framei,2,jj), 'r+', 'MarkerSize', 12, 'LineWidth',1);
                    ellipse(this.areas(framei,jj).MajorAxisLength/2, this.areas(framei,jj).MinorAxisLength/2, ...
                            this.orient(framei,jj), this.COM(framei,1,jj),this.COM(framei,2,jj),'r');
                    %text(this.COM(framei,1,jj)+10, this.COM(framei,2,jj), num2str(this.blobID(framei, jj)), 'Color','r', 'FontSize', 14);
                end
                line(this.bodyCOM(framei,1), this.bodyCOM(framei,2), 'Marker', 'o', 'Color', 'c','MarkerSize', 12, 'LineWidth',2);
                    %line(this.areas(framei).Extrema(:,1), this.areas(framei).Extrema(:,2), 'Marker', '.', 'Color', 'c');
                    %[u, v] = pol2cart(this.orient(framei), this.vel(framei)*.1);
                    %quiver(this.COM(framei,1), this.COM(framei,2), u,v, 'LineWidth', 2); %plots an orientation arrow
                if this.tailVisible(framei)
                    line(this.areas(framei,this.tailblob(framei)).Centroid(1), this.areas(framei,this.tailblob(framei)).Centroid(2), ...
                        'Marker', 'x', 'Color', 'm','MarkerSize', 12, 'LineWidth',2);
                end
                if (~isnan(this.noseblob(framei)))
                    line(this.areas(framei,this.noseblob(framei)).Centroid(1), this.areas(framei,this.noseblob(framei)).Centroid(2), ...
                        'Marker', '+', 'Color', 'g','MarkerSize', 12, 'LineWidth',2);
%                    [xv,yv] = pol2cart(this.orient(framei), this.vel(framei));
%                    quiver(this.areas(framei, this.noseblob(framei)).Centroid(1), this.areas(framei,this.noseblob(framei)).Centroid(2), xv, yv, 0);
                end
            end
           set(gca, 'xlim', [dispCrop(1) dispCrop(3)], 'ylim', [dispCrop(2) dispCrop(4)]);
        end
        %------------------------------------------------------------------
        function rgbIm = gray2rgb(this, grayIm)
            rgbIm = cat(3, grayIm, grayIm, grayIm);
        end
        
        %--------------------------------------------------------------------------------------------
        function pathIm = plotPathsOnFrame(this,inFrame,pathNums, varargin)
            if nargin < 4
                blend = 1;
            else 
                blend = varargin{1};
            end
            pathIm = inFrame;
            if isempty(this.paths)
                return;
            else
                if (~exist('pathNums','var'))  
                    pathNums = 1:length(this.paths);
                elseif(~isempty(pathNums)) %if we're given a set of paths, use that, otherwise plot all of them
                    pathNums = intersect(1:length(this.paths), pathNums); %make sure we don't try to plot anything not there
                else
                    pathNums = 1:length(this.paths);
                end
                %pathIm = cat(3, this.avgFrame, this.avgFrame, this.avgFrame);
                pathIm = cat(3, inFrame, inFrame, inFrame);
                %pathIm = zeros(this.height, this.width, 3);
                color_order = [2 1 3];
                for ii = 1:length(pathNums)
                    %set color layer to 255
                    pi = pathNums(ii);
                    c = color_order(mod(pi-1, 3)+1);
                    layer = pathIm(:,:,c);
                    orig = layer(this.paths(pi).PixelIdxList);
                    layer(this.paths(pi).PixelIdxList) = floor((orig*(1-blend) + 255*blend)./2);
                    pathIm(:,:,c) = layer;
                    oc = find(1:3 ~= c);
                    %set the other layers to 0
                    layer = pathIm(:,:,oc(1));
                    layer(this.paths(pi).PixelIdxList) = floor(orig*((2-blend)/2));
                    pathIm(:,:,oc(1)) = layer;
                    layer = pathIm(:,:,oc(2));
                    layer(this.paths(pi).PixelIdxList) = floor(orig*((2-blend)/2));
                    pathIm(:,:,oc(2)) = layer;
                end
            end
        end
        
        % ------------------------------------------------------------------------------------------------------
        %           
        function res = readFrames(this, frames, flag)
        % function vid_struct = readFrames(this, frames, flag)
        % 
        % Flag specifies if the reading is a SINGLE frame, a CONTINUOUS range, or a DISCONTINUOUS set of frames
            adjFrames = frames+double(this.frameRange(1))-1;
            if strcmp(flag, 'single')
                res = this.readerObj.read(adjFrames(1));
                res = squeeze(res(:,:,1));
            elseif strcmp(flag,'continuous')
                res = this.readerObj.read(adjFrames);
                res = squeeze(res(:,:,1,:));
            elseif strcmp(flag, 'discontinuous') %due to limitations in videoReader class, have to loop for this
                totalFrames = sum(diff(frames')'+1);
                res = zeros(this.nativeHeight, this.nativeWidth, totalFrames); %4D matrix
                res_ind=0;
                for ii = 1:size(adjFrames,1) %loop that reads
                    ind_n = adjFrames(ii,2)-adjFrames(ii,1)+1;
                    inds = (1:ind_n) + res_ind;
                    res_ind = inds(end);
                    fi = adjFrames(ii,:);
                    if ind_n == 1
                        fi = adjFrames(ii,1);
                    end
                    temp = this.readerObj.read(fi);
                    res(:,:,inds) = squeeze(temp(:,:,1,:));
                    %res(:,:,:,inds) = this.readerObj.read(fi);
                end
            end
            res = this.applyCrop(res); %crop
        end
        
        % ------------------------------------------------------------------------------------------------------
        function [movieArray, fc_val] = readMovieSection(this, frames, movieType, varargin)
        % function movieArray = readMovieSection(this, frames, movieType)
        % 
        % reads a section of movie specified by the frames list.
        % Inputs: frames - a list of frame ranges, formatted like [1 40; 80 100]
        %         movieType - string either 'orig' for original movie, 'diff' for the difference, or 'bin' for the
        %         binary
        %         varargin - a threshold value if you want to specify that for a binary movie 
            frames = unique(frames(:));
            [frame_ints, flag] = this.listToIntervals(frames);
            rawArray = this.readFrames(frame_ints, flag);
            if (size(rawArray,3) ~= length(frames))
                disp('In MouseTrackerKF.readMovieSection: readFrames has returned an array of different length than requested');
                s3 = size(rawArray,3);
                frames = 1:s3;
            end
            [movieArray, ~, ~, fc_val] = this.processFrameArray(rawArray, 1:length(frames), this.avgFrame, movieType, varargin{:});
        end
        
        % ------------------------------------------------------------------------------------------------------
        function movieArray = applyCrop(this, movieArray)
            % The crop is the number of px to cut in from the edge in this format [left, right, top, bottom]
            left = (this.crop(1)+1);
            horiz = (1:this.width) + left - 1;
            top = (this.crop(3)+1);
            vert = (1:this.height) + top - 1;
            movieArray = movieArray(vert, horiz, :,:);
            
        end
        
        % --------------------------------------------------
        function findMovieFile(this, vidFN)
            % function findMovieFile(this)
            %
            % Exists to find the movie file if the saved .mat file for the
            % tracking is moved since the VideoReader object saves the absolute
            % filename of the movie.
            
            disp(['Opening new video reader object for: ' vidFN]);
            this.readerObj = VideoReader(vidFN);
        end
        
        % ---------------------------------------------------
        function fakePath = makeMirrorPath(this,path)  
            fakePath.Area = path.Area;
            PixelList = [this.width - path.PixelList(:,1)+1, this.height - path.PixelList(:,2)+1];
            fakePath.PixelIdxList= this.height*(PixelList(:,1)-1) + PixelList(:,2);
            fakePath.PixelList = PixelList;
        end
        
        % ---------------------------------------------------
        function rethreshold(this, frameRange, p_mouse)
            % Setting the threshold differently in order to improve blob identification
        
            % There are two ways that we can specify thresholds - either directly or using the p_mouse
            % variable, which sets the threshold so that a p_mouse proportion of pixels are above it. 
            % p_mouse takes precedence.   
            
            this.p_mouse = p_mouse;
            this.findMouse(frameRange);
        end
        % ---------------------------------------------------------------------------
        function [ret_mov, avg_frame, thresh, fc_val] = processFrameArray(this, rawArray, frame_range, subFrame, movieType, varargin)
            % returns a  movie using the frames specified in the input. The movie can appear different and
            % is specified with MOVIETYPE: 'orig' gives the original movie, 'diff' provides a version with a 
            % frame subtracted, by default the mean frame, 'bin' is a binary thresholded image. There is
            % an optional argument 'subFrame' specifying a frame to subtract from each frame (leave empty, [], 
            % if you don't want to specify) in order to improve the thresholding of certain objects. If giving multiple
            % thresholds, then they will be in the 4th dimension of the array.
            % RAWARRAY is in the format of dim 1,2 - height,width, 3- frame#.
            % VARARGIN{1} is the threshold level(s) for the movie, leaving
            % it out gives the default, and the returned movie is a cell
            % area of 
           
            
            %Image Processing Settings
            %thresh(1) = .1; % the threshold level
            p_mouse = this.p_mouse;
            %p_mouse = .0007; %the prior probability of a mouse pixel.  Influences the threshold.
            
            erode_size = 3; %the size of erosion mask
            do_hpfilter = 1; %flag for highpass filtering
            alpha = .5; %The parameter for an unsharp filter - subtracts a blurred image from the image to sharpen original
            new_mov = rawArray(:,:,frame_range);
            if ~isempty(this.fcArea) % excludes the area used for counting frames
                fca = this.fcArea;
                fc_mov = new_mov(fca(2):fca(4), fca(1):fca(3), :);
                fc_val = squeeze(sum(sum(fc_mov,1),2));
                new_mov(fca(2):fca(4), fca(1):fca(3), :) = 0;
            else
                fc_val = NaN*ones(length(frame_range), 1);
            end
            % make a movie from the average frame to subtract
            if ~isempty(subFrame)
                avg_frame = subFrame;
            else
                avg_frame = uint8(round(mean(new_mov,3)));
            end
            avg_mov = uint8(repmat(avg_frame, [1 1 size(new_mov,3)]));
            %diff_mov = imabsdiff(new_mov, avg_mov); %this should give a nice moving blob.
            diff_mov = new_mov - avg_mov; %this should give a nice moving blob.
            if do_hpfilter
                h = fspecial('unsharp', alpha);
                diff_mov = imfilter(diff_mov, alpha);
            end
            if(this.boostContrast)
                diff_mov = increaseMovContrast(diff_mov);
            end
            % Set the threshold for making a binary image
            thresh = this.default_thresh; %.08-.12 have worked well after image normalization
            if ~isempty(varargin) && ~isempty(varargin{1}) % I'm not exactly sure MATLAB is making empty cells
                thresh = varargin{1}; 
                if nargin > 6 %p_mouse input
                    p_mouse = varargin{2};
                end
                if thresh == 0
                    % The way we are determining the threshold value is to take the brightest p_mouse
                    % proportion of pixels
                    sorted = sort(diff_mov(:), 'descend');
                    thresh_i = round(p_mouse*length(sorted));
                    thresh = sorted(thresh_i);
                    %thresh = .6 * graythresh(diff_mov(:)); % A way of doing auto thresholding
                    this.used_thresh = thresh;
                end
            end
            if strcmp(movieType, 'bin')
                ret_mov = [];
                for kk = 1:length(thresh)
                    bin_mov = diff_mov > thresh(kk);
                    %bin_mov = false(size(diff_mov));
                    nFrames = size(bin_mov,3);
                    for ii = 1:nFrames
                        %this is to get rid of the jagged edges due to something regarding movie compression, but also
                        %removes one pixel around the edge of the contiguous sections of blmtob
                        bw = bin_mov(:,:,ii);
                        %bw = im2bw(diff_mov(:,:,ii),thresh(kk));
                        bw = imerode(bw, strel('square',erode_size)); 
                        bin_mov(:,:,ii) = bw;
                    end
                    ret_mov = cat(4, ret_mov, bin_mov);
                end
            elseif strcmp(movieType, 'diff')
                ret_mov = diff_mov;
            else % show the original movie
                if(this.boostContrast)
                    new_mov = increaseMovContrast(new_mov);
                end
                ret_mov = new_mov;
            end
        end
        

    end %Methods
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STATIC METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Static, Access = private)
        function  mov_struct = convertToGray(mov_struct)
            % function that converts the movie, frame by frame into grayscale - if it's not already
            for ii=1:length(mov_struct.frames)
                mov_struct.frames(ii).cdata = rgb2gray(mov_struct.frames(ii).cdata);
            end
        end
        
        % ---------------------------------------------------------------------------
        function avg_frame = averageFrame(mov_struct, frame_range) 
            % function to that takes the specific movie structure and computes the average frame of the specified frames
            mov = mov_struct.frames(frame_range);
            new_mov = zeros([size(mov(1).cdata,1) size(mov(1).cdata,2) length(mov)], 'uint8');
            for ii=1:length(mov)
                new_mov(:,:,ii) = squeeze(mov(ii).cdata(:,:,1));
            end
            avg_frame = uint8(round(mean(new_mov,3)));
        end
        
        
        
        % ---------------------------------------------------------------------------
        function [intervals, desc] = listToIntervals(list)
        % function intervals = listToIntervals(list)
        %
        % This function breaks up a possibly discontinuous list 
        % into a set of bounded intervals, and a flag saying if 
        % it was continuous or not.  

            % regularize list
            list = sort(list);
            list = unique(list);

            diffs = diff(list);
            breaks = find(diffs > 1);
            intervals = zeros(length(breaks)+1, 2);
            intervals(1,1) = list(1);
            for i=1:(length(breaks))
                if ~isempty(breaks) %also will be only thing that executes in only loop
                    intervals(i,2) = list(breaks(i));
                    intervals(i+1, 1) = list(breaks(i)+1);
                end
            end
            intervals(length(breaks)+1, 2) = list(end);
            if (size(intervals) > 1)
                desc = 'discontinuous';
            else
                desc = 'continuous';
            end
        end
        
        % -------------------------------------------------
        function movieType = parseMovieType(mtype_str)
            
            if strcmpi(mtype_str, 'orig')
                movieType = 0;
            elseif strcmpi(mtype_str, 'diff')
                movieType = 1;
            elseif strcmpi(mtype_str, 'bin')
                movieType = 2;
            else
                movieType = -1;
            end
        end
              
        
    end % methods - private
end %classdef