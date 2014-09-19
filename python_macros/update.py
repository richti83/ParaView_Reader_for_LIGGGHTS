from paraview.simple import *
import glob

DIGITS = '0123456789'

def makeGlobExpression(a, b):
    #Compare two fileNames and return a glob expression that matches both

    # simple cases
    if (a == b):
        return a
    if (len(a) == 0) or (len(b) == 0):
        return '*'

    # Always make a shorter than b (saves us comparisons later on)
    if (len(a) > len(b)):
        a, b = b, a

    # Compare from the left
    left = 0
    while( (a[left] == b[left]) and (left < len(a)) ):
        left += 1

    # did we get the whole string?
    if (left == (len(a)-1)):
        return a + '*'

    # compare from the right
    right = -1
    while( (a[right] == b[right]) ):
        right -= 1
    # left == right should never happen, because we tested for equality

    # clear all digits that touch the cut
    leftStr = a[:left].rstrip(DIGITS)
    rightStr = a[right+1:].lstrip(DIGITS)

    # concatenate glob Expression
    globExpression = leftStr + '*' + rightStr
    return globExpression


def try_int(s):
        "Convert to integer if possible."
        try: return int(s)
        except: return s


def natsort_key(s):
        "Used internally to get a tuple by which s is sorted."
        import re
        return map(try_int, re.findall(r'(\d+|\D+)', s))


def natcmp(a, b):
        "Natural string comparison, case sensitive."
        return cmp(natsort_key(a), natsort_key(b))


def natcasecmp(a, b):
        "Natural string comparison, ignores case."
        return natcmp(a.lower(), b.lower())

# Iterate over all sources
endTime = 0
for description, source in GetSources().iteritems():
    print source.__class__.__name__
    if (source.__class__.__name__ not in ['Transform', 'ProgrammableFilter','Threshold','Text']):
     fileNames = source.FileNames
     # If there are multiple files, update them
     if (len(fileNames) > 1):
        # Estimate the glob expression
        globExpression = makeGlobExpression(fileNames[0], fileNames[1])
        # update the list
        files = glob.glob(globExpression) 
        files.sort(natcasecmp)
        source.FileNames=files
        endTime = max(endTime, len(files)-1)

        # If we are dealing with particles, set their representation to point sprites
        if (source.__class__.__name__ in ['liggghts_Reader', 'liggghts_binReader']):
            SetDisplayProperties(source, Representation="Points")
            DataRepresentation = Show(source)
            DataRepresentation.PointSpriteDefaultsInitialized = 1
            DataRepresentation.Texture = []
            DataRepresentation.RadiusTransferFunctionEnabled = 1
            DataRepresentation.RadiusMode = 'Scalar'
            DataRepresentation.Representation = 'Point Sprite'
            DataRepresentation.RadiusArray = [None, 'radius']
            DataRepresentation.RadiusIsProportional = 1
            SetActiveSource(source)

# Jump to end
print endTime
RenderView = GetRenderView()
AnimationScene = GetAnimationScene()
RenderView.ViewTime = endTime
AnimationScene.AnimationTime = endTime
Render()
