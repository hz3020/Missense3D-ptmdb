# from django.db import models

class Segment:
    def __init__(self, start, end, color):
        self.start = start
        self.end = end
        self.color = color

    def to_dict(self):
        return {
            'start': self.start,
            'end': self.end,
            'color': self.color,
        }

class SegmentList:
    def __init__(self):
        self.segments = []

    def add_segment(self, start, end, color):
        segment = Segment(start, end, color)
        self.segments.append(segment)

    def remove_segment(self, index):
        if 0 <= index < len(self.segments):
            del self.segments[index]
        else:
            print("Index out of range")

    def get_segments(self):
        return self.segments

    def __repr__(self):
        return "\n".join(str(segment) for segment in self.segments)

class Spagetti:
    def __init__(self, residue, position, position2, color):
        self.position = position
        self.position2 = position2
        self.residue = residue
        self.color = color

    def to_dict(self):
        return {
            'residue':self.residue,
            'position': self.position,
            'position2': self.position2,
            'color': self.color,
        }
        
class SpagettiList:
    def __init__(self):
        self.spagettis = []

    def add_spagetti(self, residue, position, position2, color):
        spagetti = Spagetti(residue, position, position2, color)
        self.spagettis.append(spagetti)


    def get_spagettis(self):
        return self.spagettis

    def __repr__(self):
        return "\n".join(str(spagetti) for spagetti in self.spagettis)

    

 

