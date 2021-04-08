# Help source:
#   * https://datacarpentry.org/image-processing/08-edge-detection/
#   * https://stackoverflow.com/questions/61432335/blob-detection-in-python
#   * https://pythonexamples.org/python-opencv-cv2-resize-image/
#   *

import cv2
import numpy as np
import pygame
import os

# Constants

scale_procent = 0.3

#
# THIS IS FOR EVERY IMAGE
#

# THIS FUNCTIONS AUTOMATICLY DOES THE WORK FOR THE USER FOR THE REST OF THE FILES
def ImageAlteration(image_location, scale_procent, zoom, lowest_point, axies_of_rotation): # ADD THINGS SUCH AS SCALEABILITY AND STUFF

    """
    image_location      The source of the image
    scale_procent       To rescale the image to a proper size
    zoom                Four points in which the ball is within => Change name to zoom_points
    lowest_point        The lowest point of the ball
    axies_of_rotation
    """

    print(image_location)

    image = cv2.imread(image_location)

    height = int(image.shape[0] * scale_procent)
    width =  int(image.shape[1] * scale_procent)
    image = cv2.resize(image, (width, height))

    # Zooms in to the image to avoid unnessesary background information
    image = image[zoom[2] : lowest_point, zoom[0]:zoom[0] + zoom[1]]

    # Resize the image so the rotational axies is in the middle of the image.

    right_side = width - axies_of_rotation
    left_side = width - right_side

    if left_side < right_side:
        image = image[0 : height, 0 : 2 * left_side]
    elif left_side > right_side:
        image = image[0 : height, width - 2 * right_side : width]

    # convert img to grayscale
    gray_image = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY) # MOVED UPWARDS FOR TESTING

    # do adaptive threshold on gray image
    thresh_image = cv2.adaptiveThreshold(gray_image, 1000, cv2.ADAPTIVE_THRESH_MEAN_C, cv2.THRESH_BINARY, 201, 3)

    # apply morphology open then close
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5,5))
    blob_image = cv2.morphologyEx(thresh_image, cv2.MORPH_OPEN, kernel)
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (9,9))
    blob_image = cv2.morphologyEx(blob_image, cv2.MORPH_CLOSE, kernel)


    # invert blob
    blob_image = (255 - blob_image)


    height = blob_image.shape[0]
    width  = blob_image.shape[1]
    for i in range(0, height):
        for j in range(0, width):
            n = blob_image[i, j]

            if n == 0:
                s = image[i, j].tolist()
                if s[0] < 140 and s[1] < 140 and s[2] < 140:
                    blob_image[i, j] = 255


    return blob_image

#
# USE AN EXAMPLE IMAGE TO GATHER INFORMATION ABOUT IMAGE MANIPULATION
#

def ExampleImage(image_location):

    def get_zoom(image):

        def opencv_to_pygame_image(image):
            return pygame.image.frombuffer(image.tostring(), image.shape[1::-1], "BGR")

        pygame.init()

        background = opencv_to_pygame_image(image)
        screen = pygame.display.set_mode(background.get_rect().size)
        x_max = background.get_rect().size[0]
        y_max = background.get_rect().size[1]

        zoom = 0
        zoom_specifics = [0, 0, 0, 0]
        running = True
        while running:
            for event in pygame.event.get():
                if event.type == pygame.QUIT:
                    done = True

                if event.type == pygame.MOUSEBUTTONDOWN:
                    running = False

                screen.blit(background, (0, 0))

                mouse_position = pygame.mouse.get_pos()

                zoom = mouse_position[0] / background.get_rect().size[0]

                # This is relative to the image where zero is at the bottom and max value at the top, respectivly min to the left and max to right on x-axies.
                xa = int(x_max / 2 * zoom)
                ya = int(y_max / 2 * zoom)

                xb = int(2 * (x_max / 2 - xa))
                yb = int(2 * (y_max / 2 - ya))

                zoom_specifics = [xa, xb, ya, yb]

                # Draw lines on screen
                pygame.draw.rect(screen, (0, 0, 0), (xa, ya, xb, yb), 3)

                pygame.display.update()


        pygame.display.quit()

        return zoom_specifics

    def low_point_selector(image):

        def opencv_to_pygame_image(image):
            return pygame.image.frombuffer(image.tostring(), image.shape[1::-1], "RGB")

        pygame.init()

        background = opencv_to_pygame_image(image)
        screen = pygame.display.set_mode(background.get_rect().size)
        x_max = background.get_rect().size[0]
        y_max = background.get_rect().size[1]

        low_point = 0
        running = True
        while running:
            for event in pygame.event.get():
                if event.type == pygame.QUIT:
                    done = True

                if event.type == pygame.MOUSEBUTTONDOWN:
                    running = False

                screen.blit(background, (0, 0))

                mouse_position = pygame.mouse.get_pos()

                low_point = mouse_position[1]
                pygame.draw.line(screen, (0, 0, 0), (0, mouse_position[1]), (x_max, mouse_position[1]))

                pygame.display.update()

        pygame.display.quit()
        return low_point

    def axies_of_rotation(image):
        def opencv_to_pygame_image(image):
            return pygame.image.frombuffer(image.tostring(), image.shape[1::-1], "RGB")

        pygame.init()

        background = opencv_to_pygame_image(image)
        screen = pygame.display.set_mode(background.get_rect().size)
        x_max = background.get_rect().size[0]
        y_max = background.get_rect().size[1]

        axies = 0
        running = True
        while running:
            for event in pygame.event.get():
                if event.type == pygame.QUIT:
                    done = True

                if event.type == pygame.MOUSEBUTTONDOWN:
                    running = False

                screen.blit(background, (0, 0))

                mouse_position = pygame.mouse.get_pos()

                axies = mouse_position[0]
                pygame.draw.line(screen, (0, 0, 0), (mouse_position[0], 0), (mouse_position[0], y_max))

                pygame.display.update()

        pygame.display.quit()
        return axies


    # read image => CHANGE THIS TO EXAMPLEFILE => SORT NAMES AND USE LOWEST ANGLE AS REFERANCE
    image = cv2.imread(image_location)

    height = int(image.shape[0] * scale_procent)
    width =  int(image.shape[1] * scale_procent)
    image = cv2.resize(image, (width, height))

    # Zooms in to the image with choosen zoomin
    zoom = get_zoom(image)
    lowest_point = low_point_selector(image)
    axies = axies_of_rotation(image)

    del image

    return zoom, lowest_point, axies

def main(source, target):

    if source[::-1][0] == '/':
        source = '/'.join(source.split('/')[:len(source.split('/')) - 1])

    """
    source = directory of all the raw images
    target = directory to save all the processed images
    """
    print("\nIMAGE PROCESSING:")

    # Get all images from the directory
    from os import listdir
    from os.path import isfile, join
    onlyfiles = [f for f in listdir(source) if isfile(join(source, f))]
    print("Folder:", source)
    print("Images:", onlyfiles)

    zoom, lowest_point, axies = ExampleImage(source + '/' + onlyfiles[0])

    for image in onlyfiles:
        print("Analyzing", image)
        analyzed_image = ImageAlteration(source + '/' + image, scale_procent, zoom, lowest_point, axies)
        cv2.imwrite(target + '/' + image, analyzed_image)


# cv2.imshow("RESULT", result_image )
# cv2.waitKey(0)
