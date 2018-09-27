/*
 * Copyright 2014 MovingBlocks
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package org.terasology.math.geom;

import com.google.common.collect.Sets;
import org.junit.Test;

import java.util.Arrays;
import java.util.List;
import java.util.Set;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

public class Region3iTest {

    @Test
    public void testCreateRegionWithMinAndSize() {
        List<Vector3i> mins = Arrays.asList(new Vector3i(), new Vector3i(1, 1, 1), new Vector3i(3, 4, 5));
        List<Vector3i> size = Arrays.asList(new Vector3i(1, 1, 1), new Vector3i(3, 3, 3), new Vector3i(8, 5, 2));
        List<Vector3i> expectedMax = Arrays.asList(new Vector3i(1, 1, 1), new Vector3i(4, 4, 4), new Vector3i(11, 9, 7));
        for (int i = 0; i < mins.size(); ++i) {
            Region3i region = new Region3i().setMinSize(mins.get(i), size.get(i));
            assertEquals(mins.get(i), region.min());
            assertEquals(size.get(i), region.size());
            assertEquals(expectedMax.get(i), region.max());
            assertFalse(region.isEmpty());
        }
    }

    @Test
    public void testCreateRegionWithMinMax() {
        Vector3i min;
        Vector3i max;
        Region3i region;

        min = new Vector3i();
        max = new Vector3i();
        region = new Region3i().setMinMax(min, max);
        assertEquals(min, region.min());
        assertEquals(max, region.max());
        assertEquals(new Vector3i(), region.size());
        assertTrue(region.isEmpty());


        min = new Vector3i(1, 1, 1);
        max = new Vector3i(3, 3, 3);
        region = new Region3i().setMinMax(min, max);
        assertEquals(min, region.min());
        assertEquals(max, region.max());
        assertEquals(new Vector3i(2, 2, 2), region.size());
        assertFalse(region.isEmpty());

        min = new Vector3i(3, 4, 5);
        max = new Vector3i(10, 8, 6);
        region = new Region3i().setMinMax(min, max);
        assertEquals(min, region.min());
        assertEquals(max, region.max());
        assertEquals(new Vector3i(7, 4, 1), region.size());
        assertFalse(region.isEmpty());

    }

    @Test
    public void testCreateRegionWithBounds() {
        Region3i expectedRegion = new Region3i().setMinMax(new Vector3i(-2, 4, -16), new Vector3i(4, 107, 0));

        assertEquals(expectedRegion, Region3i.createBounded(new Vector3i(-2, 4, -16), new Vector3i(4, 107, 0)));
        assertEquals(expectedRegion, Region3i.createBounded(new Vector3i(4, 4, -16), new Vector3i(-2, 107, 0)));
        assertEquals(expectedRegion, Region3i.createBounded(new Vector3i(-2, 107, -16), new Vector3i(4, 4, 0)));
        assertEquals(expectedRegion, Region3i.createBounded(new Vector3i(-2, 4, 0), new Vector3i(4, 107, -16)));
        assertEquals(expectedRegion, Region3i.createBounded(new Vector3i(4, 107, -16), new Vector3i(-2, 4, 0)));
        assertEquals(expectedRegion, Region3i.createBounded(new Vector3i(4, 4, 0), new Vector3i(-2, 107, -16)));
        assertEquals(expectedRegion, Region3i.createBounded(new Vector3i(-2, 107, 0), new Vector3i(4, 4, -16)));
        assertEquals(expectedRegion, Region3i.createBounded(new Vector3i(4, 107, 0), new Vector3i(-2, 4, -16)));

    }

    @Test
    public void testRegionEmptyIfMaxLessThanMin() {
        Region3i region = new Region3i().setMinMax(new Vector3i(0, 0, 0), new Vector3i(-1, 0, 0));
        assertTrue(region.isEmpty());
    }

    @Test
    public void testRegionEmptyIfSizeZeroOrLess() {
        Region3i region = new Region3i().setMinSize(new Vector3i(1, 1, 1), new Vector3i(0, 1, 1));
        assertTrue(region.isEmpty());
        region = new Region3i().setMinSize(new Vector3i(1, 1, 1), new Vector3i(1, -1, 1));
        assertTrue(region.isEmpty());
    }

    @Test
    public void testIterateRegion() {
        Vector3i min = new Vector3i(2, 5, 7);
        Vector3i max = new Vector3i(10, 11, 12);
        Region3i region = new Region3i().setMinMax(min, max);

        Set<Vector3i> expected = Sets.newHashSet();
        for (int x = min.x; x < max.x; ++x) {
            for (int y = min.y; y < max.y; ++y) {
                for (int z = min.z; z < max.z; ++z) {
                    expected.add(new Vector3i(x, y, z));
                }
            }
        }

        for (Vector3i pos : region) {
            assertTrue(expected.contains(pos));
            expected.remove(pos);
        }

        assertEquals("All vectors provided", 0, expected.size());
    }

    @Test
    public void testSimpleIntersect() {
        Region3i region1 = new Region3i().setMinMax(new Vector3i(), new Vector3i(32, 32, 32));
        Region3i region2 = new Region3i().setMinMax(new Vector3i(1, 1, 1), new Vector3i(17, 17, 17));
        assertEquals(region2, region1.intersect(region2));
    }

    @Test
    public void testNonTouchingIntersect() {
        Region3i region1 = new Region3i().setMinMax(new Vector3i(), new Vector3i(32, 32, 32));
        Region3i region2 = new Region3i().setMinMax(new Vector3i(103, 103, 103), new Vector3i(170, 170, 170));
        assertEquals(Region3i.EMPTY, region1.intersect(region2));
    }

    @Test
    public void testEncompasses() {
        Region3i region = new Region3i().setMinSize(new Vector3i(), new Vector3i(1, 1, 1));
        assertTrue(region.contains(0, 0, 0));

        assertFalse(region.contains(1, 0, 0));
        assertFalse(region.contains(1, 0, 1));
        assertFalse(region.contains(0, 0, 1));
        assertFalse(region.contains(-1, 0, -1));
        assertFalse(region.contains(-1, 0, 0));
        assertFalse(region.contains(-1, 0, -1));
        assertFalse(region.contains(0, 0, -1));

        assertFalse(region.contains(1, 1, 0));
        assertFalse(region.contains(1, 1, 1));
        assertFalse(region.contains(0, 1, 1));
        assertFalse(region.contains(-1, 1, -1));
        assertFalse(region.contains(-1, 1, 0));
        assertFalse(region.contains(-1, 1, -1));
        assertFalse(region.contains(0, 1, -1));

        assertFalse(region.contains(1, -1, 0));
        assertFalse(region.contains(1, -1, 1));
        assertFalse(region.contains(0, -1, 1));
        assertFalse(region.contains(-1, -1, -1));
        assertFalse(region.contains(-1, -1, 0));
        assertFalse(region.contains(-1, -1, -1));
        assertFalse(region.contains(0, -1, -1));
    }

    @Test
    public void testNearestPointToWhenEncompasses() {
        Region3i region = new Region3i().setMinMax(new Vector3i(), new Vector3i(4, 4, 4));
        assertEquals(new Vector3i(2, 1, 1), region.getNearestPointTo(new Vector3i(2, 1, 1)));
    }

    @Test
    public void testNearestPointToAlongSide() {
        Region3i region = new Region3i().setMinMax(new Vector3i(), new Vector3i(4, 4, 4));
        assertEquals(new Vector3i(4, 2, 1), region.getNearestPointTo(new Vector3i(15, 2, 1)));
    }

    @Test
    public void testNearestPointToAwayFromCorner() {
        Region3i region = new Region3i().setMinMax(new Vector3i(), new Vector3i(4, 4, 4));
        assertEquals(new Vector3i(4, 4, 4), region.getNearestPointTo(new Vector3i(15, 12, 7)));
    }
}