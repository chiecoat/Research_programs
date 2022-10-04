import numpy as np
import cv2 as cv
import matplotlib.pyplot as plt
import time

class Session:
    "This class encompasses the full mask creation session"

    #Gonna ask you to feed in an image directy to session,
    #However you load in the image to python is your problem
    def __init__(self, filepath,change_num = 0, act_part = 1):

    
        image = cv.imread(filepath)
        image = cv.cvtColor(image, cv.COLOR_BGR2GRAY)
        #Told you I wouldn't write exception handling,
        #but really need to put something in for when
        #filepath isn't an image (or a string)

        #Reference for display
        self.image = image

        #Make the mask to be edited
        self.mask = np.zeros(image.shape)

        #Dimensions of the mask will be useful
        self.mask_size = image.shape

        #Which particle are we currently marking
        self.act_part = act_part

        #Dictionary of changes made to mask for reference
        #Format is key = change number
        #value is list of lists, [[[a,b],[p0,p1]],...]
        #first element is list of pixel index, [y,x]
        #second element is list of pixel value change from old to new [p0,p1]
        #array format to accomodate changing multiple pixels at a time
        self.change_list = {}
        self.change_num = 0


        #need a drawing boolean to control the later drawing features
        self.drawing = False

        #need a stored value to capture drawn lines as single instances
        self.draw_start_ind = 0

        #Brush size for free hand drawing
        self.brush_size = 2

        #Alpha value for mask transparency
        self.alpha=.7

        #Variable for which drawing mode we are in
        #drawing = 0, eraser = 1, zoom = 2
        self.mode = 0

        #Copy of previous mask for memory purposes
        self.prev_mask = np.copy(self.mask)

        #Store zoom values always, set to full canvas size initially
        self.startx = 0
        self.starty = 0
        self.endx = image.shape[1]-1
        self.endy = image.shape[0]-1

        #Set the amount the window zooms when zooming
        self.zoomval = 100

        #Check if the image changed before drawing
        self.img_change = True

    #Move on to the next particle to be identified
    def IterParticle(self,void):
        last_part = np.amax(self.mask)
        self.act_part = last_part+1
        self.img_change = True

    #Undo changes made to the mask (We will check if this is a valid elsewhere)
    def Undo(self):
        prev_changes = self.change_list[self.change_num - 1] #Find the change log for the previous state


        #Check if there were multiple changes in the change log
        if type(prev_changes[0][0]) == int:
            self.mask[prev_changes[0][0],prev_changes[0][1]] = prev_changes[1][0] #Change the pixel to the previous state
        else:
            for change in prev_changes:
                self.mask[change[0][0],change[0][1]] = change[1][0] #Change the pixel to the previous state
        print(self.change_num)
        self.change_num-=1 #Set our current state backwards
        print(self.change_num)

    #Redo changes we didn't want to undo. 
    def Redo(self):
        
        future_changes = self.change_list[self.change_num] #Find the change log for the now future state

        #Check if there were multiple changes in the change log
        if type(future_changes[0][0]) == int:
            self.mask[future_changes[0][0],future_changes[0][1]] = future_changes[1][0] #Change the pixel to the previous state
        else:
            for change in future_changes:
                self.mask[change[0][0],change[0][1]] = change[1][1] #Change the pixel to the future state
        self.change_num+=1 #Set our current state to be 1 advanced, again

    #Check if there is data to redo
    def CheckRedo(self):
        if np.amax(list(self.change_list.keys()))>=self.change_num:
            return 1
        else:
            return 0
    #Check if there memory left to undo
    def CheckUndo(self):
        if self.change_num == 0:
            return 0
        elif np.amin(list(self.change_list.keys()))<self.change_num:
            return 1
        else:
            return 0


    #If a new action is taken after redo'ing, delete the forward stack of redo
    def PruneRedo(self):
        mem_steps = list(self.change_list.keys())
        print(mem_steps)

        if self.change_num < np.amax(mem_steps):
            for key in range(self.change_num,np.amax(mem_steps)+1):
                del self.change_list[key]

    #Remove stored Undo's if they are more than 50 steps behind
    def PruneUndo(self):
        mem_steps = list(self.change_list.keys())

        if len(mem_steps)>50:
            del self.change_list[np.amin(mem_steps)]


    #If you want to start from scratch
    def ClearMask(self):
        self.mask = np.zeros(self.image.shape)

    #Handle All of our drawing functionality
    def draw_circle(self,event,x,y,flags,param):
        y = y+self.starty
        x = x+self.startx
        self.img_change = True

        #Any zoom shenanigans are in here
        if self.mode == 2:
            if event == cv.EVENT_LBUTTONDOWN:
                xgap = self.endx - self.startx
                ygap = self.endy - self.starty
                zoomval = self.zoomval

                if xgap >zoomval:
                    zoom_split = ((xgap-x+self.startx)/xgap) #Calculate displacement vector from current boundary to new center
                    self.startx += int((1-zoom_split)*zoomval)
                    self.endx -= int((zoom_split)*zoomval)

                if ygap >zoomval:
                    zoom_split = ((ygap-y+self.starty)/ygap) #Calculate displacement vector from current boundary to new center
                    self.starty += int((1-zoom_split)*zoomval)
                    self.endy -= int((zoom_split)*zoomval)
                xgap = self.endx - self.startx
                ygap = self.endy - self.starty

                
            elif event == cv.EVENT_RBUTTONDOWN:
                zoomval = self.zoomval
                #Add zoom value to current boundary
                xgap = self.endx - self.startx + zoomval
                ygap = self.endy - self.starty + zoomval

                xdiff = int(xgap/2)
                ydiff = int(ygap/2)

                #Want better way to check against border cases
                #left x boundary
                if x-xdiff<0:
                    self.startx = 0
                else:
                    self.startx = x-xdiff
                #right x boundary
                if x+xdiff>self.mask_size[0]:
                    self.endx = self.mask_size[0]-1
                else:
                    self.endx = x+xdiff

                #left y boundary
                if y-ydiff<0:
                    self.starty = 0
                else:
                    self.starty = y-ydiff
                #right y boundary
                if y+ydiff>self.mask_size[1]:
                    self.endy = self.mask_size[1]-1
                else:
                    self.endy = y+ydiff
                

                
                

                
            

        elif event == cv.EVENT_LBUTTONDOWN:
            self.drawing = True
            self.draw_start_ind = self.change_num
            if self.mask[y,x] != 0:
                self.act_part = self.mask[y,x]

            #Check memory of undo/redo before making new actions
            #elif self.change_num != 0:
            #    print(self.change_list.keys())
            #    self.PruneRedo()
            #    print(self.change_list.keys())
            #    self.PruneUndo()


    
        #Start by swapping around on a given pixel
        elif event == cv.EVENT_LBUTTONUP and self.draw_start_ind == self.change_num:
            self.drawing = False
            if self.mask[y,x] == 0:
                self.mask[y,x] = self.act_part
                self.change_list[self.change_num] = [[y,x],[0,self.act_part]]
                self.change_num+=1
            elif self.mask[y,x] == self.act_part:
                self.mask[y,x] = 0
                self.change_list[self.change_num] = [[y,x],[self.act_part,0]]
                self.change_num +=1
            self.prev_mask = np.copy(self.mask)

        elif event == cv.EVENT_LBUTTONUP:
            self.drawing = False
            for memkey in range(self.draw_start_ind+1,self.change_num):
                self.change_list[self.draw_start_ind].extend(self.change_list[memkey])
                del self.change_list[memkey]
            self.change_num = self.draw_start_ind+1
            self.prev_mask = np.copy(self.mask)

            


        #Can also try and drag and paint (will add slider for brush size later). Currently in a square pattern
        elif event == cv.EVENT_MOUSEMOVE:
            if self.drawing == True:
                size = self.brush_size #number of pixels around center pixel that will be drawn, can be set as function of slider later

                #Check if we are drawing or erasing
                if self.mode == 0:
                    draw_val = self.act_part
                else:
                    draw_val = 0

                self.mask[y-size:y+size,x-size:x+size] = draw_val*np.ones((2*size,2*size))

                #Really clunky for loop to iterate through area that is changed
                changelog = []
                for yval in range(y-size,y+size):
                    changelog.extend([[[yval,xval],[self.prev_mask[yval,xval],draw_val]] for xval in range(x-size,x+size)])

                self.change_list[self.change_num] = changelog
                self.change_num+=1

        
        
                        

            
    def DrawMap(self):
        #Identify indices
        cur_pixels = self.mask == self.act_part
        old_pixels1 = self.mask >0
        old_pixels2 = self.mask!= self.act_part
        old_pixels = old_pixels1*old_pixels2

        colormask = np.dstack([self.image,self.image,self.image])
        blankmask = np.copy(colormask)

        colormask[cur_pixels] = [0,0,255]
        colormask[old_pixels] = [255,0,0]
        

        return colormask,blankmask


    

    def Memory(self,state):
        
        #When the Memory slider is changed, do appropriate action
        if state == 2:
            print('Redid')
            if self.CheckRedo(): #Check if there is a state to 'redo'
                print('no really')
                self.Redo()
        elif state == 0: #Check if there is a state to 'undo'
            if self.CheckUndo():
                self.Undo()

    #Make Zoomed in Window
    def ZoomWind(self,canvas):
        return canvas[self.starty:self.endy,self.startx:self.endx]
        

    def ChangeBrushSize(self,state):

        #When the Brush Size slider changes, update the brush size
        self.brush_size = state

    def ChangeAlpha(self,state):
        self.alpha = state/100.
        self.img_change = True

    def ChangeZoom(self,state):
        self.zoomval = state
        
    def StartDraw(self):

        #Create instance of window, catch mouse input, and display image
        cv.namedWindow("Mask Maker",cv.WINDOW_NORMAL)
        cv.resizeWindow("Mask Maker",900,800)
        
        cv.setMouseCallback("Mask Maker", self.draw_circle) 
        cv.imshow("Mask Maker", self.image) 

        #Make window for buttons/sliders
        cv.namedWindow('slider',cv.WINDOW_NORMAL)
        
        #Make 'Button' for iterating particles
        switch = 'New Particle'
        cv.createTrackbar(switch,'slider',0,1,self.IterParticle)

        #Make Undo/Redo toggle
        UnRe = '0: Undo \n1: Redo'
        cv.createTrackbar(UnRe,'slider',0,2,self.Memory)

        #Make Brush Slider
        cv.createTrackbar('Brush Size','slider',1,20,self.ChangeBrushSize)
        cv.setTrackbarPos('Brush Size','slider',self.brush_size)

        #Make Alpha Slider
        cv.createTrackbar('Opacity','slider',0,100,self.ChangeAlpha)
        cv.setTrackbarPos('Opacity','slider',int(100*self.alpha))

        #Make Zoom Slider
        cv.createTrackbar('Zoom','slider',0,300,self.ChangeZoom)
        cv.setTrackbarPos('Zoom','slider',self.zoomval)
        while True:
            if self.img_change:
                cmask,bmask = self.DrawMap()
                alpha_img = cv.addWeighted(cmask,self.alpha,bmask,1-self.alpha,0,bmask)

                #Zoom mask and image, if needed
                zoom_img = self.ZoomWind(self.image)
                zoom_mask = self.ZoomWind(bmask)
                cv.imshow("Mask Maker", zoom_img)
                cv.imshow('Mask Maker', zoom_mask)
                self.img_change = False

            #visually show the new particle has been selected
            if cv.getTrackbarPos(switch,'slider') == 1:
                time.sleep(.4)
                cv.setTrackbarPos(switch,'slider',0)

            #Visually show the undo/redo process
            if cv.getTrackbarPos(UnRe,'slider') != 1:
                time.sleep(.15)
                cv.setTrackbarPos(UnRe,'slider',1)

            #Swap Mode
            key = cv.waitKey(1) & 0xFF

            #If you click the keyboard
            if key != 255:
                self.img_change = True

            #Toggle Eraser Function
            #101 == e
            if key == 101:
                if self.mode == 1:
                    self.mode = 0
                    print('Drawing Tool is on')
                else:
                    self.mode = 1
                    print('Eraser Tool is on')

            #Toggle Zoom Function
            #122 == z
            elif key == 122:
                if self.mode == 2:
                    self.mode = 0
                    print('Drawing tool is on')
                else:
                    self.mode = 2
                    print('Zoom tool is on')

                    
                    

            #Pan while zoomed in
            #100 == d
            elif key == 100 and self.endx < self.mask_size[1]-5:
                self.startx+=5
                self.endx +=5
            #97 == a
            elif key == 97 and self.startx > 5:
                self.startx-=5
                self.endx -=5
            #115 == s
            elif key == 115 and self.endy < self.mask_size[0]-5:
                self.starty+=5
                self.endy +=5
            #119 == w
            elif key == 119 and self.starty > 5:
                self.starty-=5
                self.endy -=5
              
            elif key == 27: 
                break
                  
        cv.destroyAllWindows()

    def SaveMask(self,filepath):
        np.savetxt(filepath,self.mask)

